#include <RcppArmadillo.h>

#include "Palantir_Core/Phylogeny.hpp"
#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/GeneticCode.hpp"
#include "Palantir_Core/SiteSimulation.hpp"
#include "Palantir_Core/Simulate.hpp"
#include "Palantir_Core/Util.hpp"

#include "RcppPalantir.hpp"

#include <sys/stat.h>

using namespace Rcpp;
using namespace std;

// FIX (2026-08-20, M5): the `type` field distinguishes a branch-length tree
// ("phylogeny") from a guide tree whose branch "lengths" are mode indices
// ("mode"). It was never validated on the way in and never checked on the way
// out, so passing the branch-length tree where a guide tree was expected ran
// every branch under models[[1]] without complaint.
static void check_phylogeny_type(List phylogeny, std::string expected, std::string argument)
{
    if(!has_class(phylogeny, "Phylogeny")) {
        stop("Argument `" + argument + "` should be of class `Phylogeny`");
    }
    string type = get_attr(phylogeny, "type");
    if(type != expected) {
        stop("Argument `" + argument + "` should be a phylogeny of type `" +
             expected + "`, not `" + type + "`. Build it with "
             "Phylogeny(path, type = \"" + expected + "\").");
    }
}

// FIX (2026-08-20, M4): Phylogeny() used to accept anything ifstream could
// open -- a directory, an empty file, arbitrary text, or a file holding several
// trees (only the first was used). The newick parser then either built a
// nonsense tree or failed with an unrelated message. These are the cheap
// structural checks; they are deliberately not a full grammar.
static void validate_newick(const std::string& path, const std::string& newick)
{
    string trimmed(newick);
    removeWhitespace(trimmed);

    if(trimmed.empty()) {
        stop("File '" + path + "' is empty; expected a newick tree");
    }

    long depth = 0;
    for(const char& c : trimmed) {
        if(c == '(') {
            depth++;
        } else if(c == ')') {
            depth--;
            if(depth < 0) {
                stop("File '" + path + "' is not a valid newick tree: "
                     "unbalanced parentheses (a ')' with no matching '(')");
            }
        }
    }
    if(depth != 0) {
        stop("File '" + path + "' is not a valid newick tree: unbalanced "
             "parentheses (" + std::to_string(depth) + " unclosed '(')");
    }

    size_t terminator = trimmed.find(';');
    if(terminator == string::npos) {
        stop("File '" + path + "' is not a valid newick tree: it does not end "
             "with ';'");
    }
    if(terminator != trimmed.size() - 1) {
        stop("File '" + path + "' contains more than one tree; PalantiR reads "
             "a single newick tree per file");
    }
}

// [[Rcpp::export]]
List Phylogeny(std::string newick_path, std::string type = "phylogeny")
{
    if(type != "phylogeny" && type != "mode") {
        stop("Argument `type` should be either \"phylogeny\" (branch lengths) "
             "or \"mode\" (a guide tree whose branch lengths are mode "
             "indices), not \"" + type + "\"");
    }

    struct stat info;
    if(stat(newick_path.c_str(), &info) != 0) {
        stop("Could not find file " + newick_path);
    }
    if(S_ISDIR(info.st_mode)) {
        stop("'" + newick_path + "' is a directory, not a newick file");
    }

    string newick;
    ifstream f(newick_path);
    if(f.good()) {
        newick = file_to_string(newick_path);
    } else {
        stop("Could not read file " + newick_path);
    }
    validate_newick(newick_path, newick);

    Palantir::Phylogeny tree(newick);

    List phylo = List::create(
        _["json"] = tree.to_JSON(),
        _["string"] = tree.to_string(),
        _["newick"] = newick,
        _["n_nodes"] = tree.n_nodes,
        _["type"] = type
    );
    phylo.attr("class") = "Phylogeny";
    return phylo;
}

// FIX (2026-08-20, C7): the engine RNG is a single mt19937 that is independent
// of R's generator, so R's set.seed() has no effect on it. Call this instead to
// make a run of PalantiR simulations reproducible; unseeded sessions keep
// drawing their seed from std::random_device.
// [[Rcpp::export]]
void set_palantir_seed(double seed)
{
    unsigned long long s = as_count(seed, "seed", 0);
    if(s > 4294967295ULL) {
        stop("Argument `seed` should be between 0 and 4294967295");
    }
    Palantir::seed_rng((unsigned int) s);
}

// FIX (2026-08-20, M7): the valid genetic code names were only discoverable by
// reading the C++ table; use_genetic_code() accepted anything and the failure
// surfaced much later. The R layer validates against this list.
// [[Rcpp::export]]
std::vector<std::string> genetic_code_names()
{
    std::vector<std::string> names;
    for(const auto& code : Palantir::GeneticCode::GeneticCode_table) {
        names.push_back(code.first);
    }
    return names;
}

// FIX (2026-08-20, M2): `nucleotide_equilibrium` is optional and defaults to
// NULL, which reproduces the previous (degeneracy-tilted) inversion exactly.
// When it is supplied, the mutational weight of each amino acid's codons is
// divided out so that a MutationSelection model built on the returned fitness
// vector -- with the SAME nucleotide equilibrium -- has the requested
// amino-acid marginal.
// [[Rcpp::export]]
arma::vec equilibrium_to_fitness(
    arma::vec equilibrium,
    double population_size,
    Rcpp::Nullable<Rcpp::NumericVector> nucleotide_equilibrium = R_NilValue)
{
    // FIX (2026-08-20): population_size here is a real effective size (the
    // formula is smooth in N and callers pass branch-length-weighted means);
    // require positive and finite, not integral.
    if (!std::isfinite(population_size) || population_size <= 0) {
        stop("Argument `population_size` should be a positive finite number, not "
             + std::to_string(population_size));
    }
    double N = population_size;

    if(nucleotide_equilibrium.isNull()) {
        return Palantir::equilibrium_to_fitness(equilibrium, N);
    }

    arma::vec nucleotide = as<arma::vec>(NumericVector(nucleotide_equilibrium));
    if(nucleotide.n_elem != Palantir::Nucleotide::size) {
        stop("Argument `nucleotide_equilibrium` should have " +
             std::to_string(Palantir::Nucleotide::size) + " elements (T, C, A, "
             "G), not " + std::to_string(nucleotide.n_elem));
    }
    if(!nucleotide.is_finite() || nucleotide.min() <= 0) {
        stop("Argument `nucleotide_equilibrium` should be strictly positive");
    }
    check_length(equilibrium.n_elem, Palantir::AminoAcid::size, "equilibrium");

    Palantir::GeneticCode g(get_genetic_code_name());
    return Palantir::equilibrium_to_fitness(equilibrium, N, nucleotide, g);
}

// [[Rcpp::export]]
List simulate_over_phylogeny(
    List phylogeny,
    List model,
    List sequence,
    double rate = 1)
{
    // Type checking
    check_phylogeny_type(phylogeny, "phylogeny", "tree");
    if(!has_class(model, "SubstitutionModel")) {
        stop("Argument `substitution_model` should be of class `SubstitutionModel`");
    }
    if(!has_class(sequence, "Sequence")) {
        stop("Argument `sequence` should be of class `Sequence`");
    }

    arma::mat transition = model["transition"];
    arma::mat sampling = model["sampling"];

    string newick = phylogeny["newick"];
    Palantir::Phylogeny p(newick);
    uvec states = sequence["index"];

    vector<Palantir::SiteSimulation> sims =
        Palantir::Simulate::sequence_over_phylogeny(
            p, transition, sampling, states, rate);

    List substitutions = site_simulations_to_list(sims, p);
    decorate_substitutions(substitutions, model["type"]);

    CharacterMatrix alignment = get_alignment(sims, p, model["type"]);

    List simulation = List::create(
        _["phylogeny"] = phylogeny,
        _["model"] = model,
        _["substitutions"] = DataFrame::create(substitutions, _["stringsAsFactors"] = false),
        _["alignment"] = alignment,
        // FIX (2026-08-20, MIN6): `NULL` here is the C++ null pointer literal,
        // which Rcpp wraps as the integer 0, so sim$intervals was numeric(0-ish)
        // and is.null(sim$intervals) was FALSE for every homogeneous
        // simulation. R_NilValue is R's NULL.
        _["intervals"] = R_NilValue,
        _["type"] = model["type"]
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}

//[[Rcpp::export]]
DataFrame phylogeny_to_intervals(List phylogeny, List mode_phylogeny)
{
    // FIX (2026-08-20, M5/MIN7): validate before dereferencing, and check the
    // phylogeny types rather than trusting the caller to pass the guide tree in
    // the second slot.
    check_phylogeny_type(phylogeny, "phylogeny", "phylogeny");
    check_phylogeny_type(mode_phylogeny, "mode", "mode_phylogeny");

    string newick = phylogeny["newick"];
    Palantir::Phylogeny p(newick);

    string mode_newick = mode_phylogeny["newick"];
    Palantir::Phylogeny pa(mode_newick);

    vector<Palantir::IntervalHistory> tree_intervals = p.to_intervals(pa);
    List intervals = interval_histories_to_list(tree_intervals, p);

    return DataFrame(intervals);
}

//[[Rcpp::export]]
List simulate_over_interval_phylogeny(
    List phylogeny,
    List mode_phylogeny,
    List models,
    List sequence,
    unsigned long long start_mode,
    double rate = 1,
    double segment_length = 0.001,
    double tolerance = 0.001,
    std::string rescale_method = "segments")
{

    // Type checking
    // FIX (2026-08-20, MIN7): models[0] used to be dereferenced for its `type`
    // BEFORE anything was validated, so an empty list or a list of non-models
    // failed with Rcpp's "index out of bounds" / "Object was created without
    // names" instead of a message naming the argument.
    // FIX (2026-08-20, M5): the guide tree's `type` is now required to be
    // "mode"; passing the branch-length tree used to run every branch under
    // models[[1]] silently.
    check_phylogeny_type(phylogeny, "phylogeny", "tree");
    check_phylogeny_type(mode_phylogeny, "mode", "mode_tree");
    if(models.size() == 0) {
        stop("Argument `substitution_models` should contain at least one model");
    }
    for(ullong i = 0; i < models.size(); i++) {
        List s = models[i];
        if(!has_class(s, "SubstitutionModel")) {
            stop("Each argument in `substitution_models` should be of class `SubstitutionModel`");
        }
    }

    List first_model = models[0];
    string model_type = get_attr(first_model, "type");

    for(ullong i = 0; i < models.size(); i++) {
        List s = models[i];
        if(get_attr(s, "type") != model_type) {
            stop("All models in `substitution_models` should have the same `type`");
        }
        if(get_attr(s, "scaling_type") != get_attr(first_model, "scaling_type")) {
            stop("All models in `substitution_models` should have the same `scaling_type`");
        }
    }
    if(!has_class(sequence, "Sequence")) {
        stop("Argument `sequence` should be of class `Sequence`");
    }

    Palantir::GeneticCode g(get_genetic_code_name());

    string newick = phylogeny["newick"];
    Palantir::Phylogeny p(newick);

    string mode_newick = mode_phylogeny["newick"];
    Palantir::Phylogeny pa(mode_newick);

    vector<Palantir::IntervalHistory> tree_intervals = p.to_intervals(pa);
    List intervals = interval_histories_to_list(tree_intervals, p);

    if(!compare_modes(models, intervals)) {
        stop("The modes on the `mode_tree` phylogeny should index the elements in `substitution_models` list");
    }

    uvec states = sequence["index"];
    uvec codons = Palantir::CompoundCodon::to_codon_index(states, g);

    vector<vec> equilibrium;
    vector<mat> transition;
    vector<mat> sampling;

    for(ullong i = 0; i < models.size(); i++) {
        List substitution_model = models[i];
        equilibrium.push_back(substitution_model["equilibrium"]);
        transition.push_back(substitution_model["transition"]);
        sampling.push_back(substitution_model["sampling"]);
    }

    // FIX (2026-08-17): scaling_type was never passed, so the segment rescaler
    // always used its default of "synonymous" whatever the models were built
    // with (PROJECT-RECORD 8E.5a). NOTE: with synonymous models this is a no-op;
    // with "substitution" models it changes behaviour, and pinning the TOTAL
    // rate during the transient is not what the analysis wants. Use synonymous.
    string scaling_type = get_attr(first_model, "scaling_type");
    vector<Palantir::SiteSimulation> sims = Palantir::Simulate::sequence_over_intervals(
        p, tree_intervals, equilibrium, transition, sampling, codons, start_mode, g, rate, segment_length, tolerance, scaling_type, rescale_method);

    List substitutions = site_simulations_to_list(sims, p);

    decorate_substitutions(substitutions, model_type);

    CharacterMatrix alignment = get_alignment(sims, p, model_type);

    List simulation = List::create(
        _["phylogeny"] = phylogeny,
        _["models"] = models,
        _["substitutions"] = DataFrame::create(substitutions, _["stringsAsFactors"] = false),
        _["alignment"] = alignment,
        _["intervals"] = DataFrame(intervals),
        _["type"] = "compound_codon"
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}

//[[Rcpp::export]]
List simulate_with_shared_substitution_heterogeneity(
    List phylogeny,
    List switching_model,
    List substitution_models,
    List sequence,
    unsigned long long start_mode,
    double rate = 1,
    double switching_rate = 1,
    double segment_length = 0.0001,
    double tolerance = 0.001)
{
    // Type checking
    // FIX (2026-08-20, M5/MIN7): validate the tree (including its `type`) and
    // every model before dereferencing substitution_models[0]. See the notes on
    // simulate_over_interval_phylogeny.
    check_phylogeny_type(phylogeny, "phylogeny", "tree");
    if(!has_class(switching_model, "SubstitutionModel")) {
        stop("Argument `switching_model` should be of class `SubstitutionModel`");
    }
    if(substitution_models.size() == 0) {
        stop("Argument `substitution_models` should contain at least one model");
    }
    for(ullong i = 0; i < substitution_models.size(); i++) {
        List s = substitution_models[i];
        if(!has_class(s, "SubstitutionModel")) {
            stop("Each argument in `substitution_models` should be of class `SubstitutionModel`");
        }
    }

    List first_model = substitution_models[0];
    string model_type = get_attr(first_model, "type");

    for(ullong i = 0; i < substitution_models.size(); i++) {
        List s = substitution_models[i];
        if(get_attr(s, "type") != model_type) {
            stop("All models in `substitution_models` should have the same `type`");
        }
        if(get_attr(s, "scaling_type") != get_attr(first_model, "scaling_type")) {
            stop("All models in `substitution_models` should have the same `scaling_type`");
        }
    }
    if(!has_class(sequence, "Sequence")) {
        stop("Argument `sequence` should be of class `Sequence`");
    }
    if(get_attr(sequence, "type") != "codon") {
        stop("Argument `sequence` should be of type `codon`");
    }

    Palantir::GeneticCode g(get_genetic_code_name());

    string newick = phylogeny["newick"];
    Palantir::Phylogeny p(newick);

    mat switching_transition = switching_model["transition"];
    mat switching_sampling = switching_model["sampling"];

    uvec states = sequence["index"];
    uvec codons = Palantir::CompoundCodon::to_codon_index(states, g);

    // One interval for each branch
    vector<Palantir::IntervalHistory> tree_intervals = Palantir::Simulate::switching_intervals(
        p, switching_transition, switching_sampling, start_mode, switching_rate);
    List intervals = interval_histories_to_list(tree_intervals, p);

    if(!compare_modes(substitution_models, intervals)) {
        stop("The modes of the `switching_model` should index the elements in `substitution_models` list");
    }

    vector<vec> equilibrium;
    vector<mat> transition;
    vector<mat> sampling;

    for(ullong i = 0; i < substitution_models.size(); i++) {
        List substitution_model = substitution_models[i];
        equilibrium.push_back(substitution_model["equilibrium"]);
        transition.push_back(substitution_model["transition"]);
        sampling.push_back(substitution_model["sampling"]);
    }

    // FIX (2026-08-17): see the note on the same call in
    // simulate_over_interval_phylogeny -- scaling_type was never passed through.
    string mm_scaling_type = get_attr(List(substitution_models[0]), "scaling_type");
    vector<Palantir::SiteSimulation> sims = Palantir::Simulate::sequence_over_intervals(
        p, tree_intervals, equilibrium, transition, sampling, codons, start_mode, g, rate, segment_length, tolerance, mm_scaling_type);

    List substitutions = site_simulations_to_list(sims, p);

    decorate_substitutions(substitutions, model_type);

    CharacterMatrix alignment = get_alignment(sims, p, model_type);

    List simulation = List::create(
        _["phylogeny"] = phylogeny,
        _["models"] = substitution_models,
        _["substitutions"] = DataFrame::create(substitutions, _["stringsAsFactors"] = false),
        _["alignment"] = alignment,
        _["intervals"] = DataFrame(intervals),
        _["type"] = "compound_codon"
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}

//[[Rcpp::export]]
List simulate_with_shared_time_heterogeneity(
    List phylogeny,
    List switching_model,
    List substitution_models,
    List sequence,
    unsigned long long start_mode,
    double rate = 1,
    double switching_rate = 1,
    double segment_length = 0.0001,
    double tolerance = 0.001)
{
    // Type checking
    // FIX (2026-08-20, M5/MIN7): as in
    // simulate_with_shared_substitution_heterogeneity above.
    check_phylogeny_type(phylogeny, "phylogeny", "tree");
    if(!has_class(switching_model, "SubstitutionModel")) {
        stop("Argument `switching_model` should be of class `SubstitutionModel`");
    }
    if(substitution_models.size() == 0) {
        stop("Argument `substitution_models` should contain at least one model");
    }
    for(ullong i = 0; i < substitution_models.size(); i++) {
        List s = substitution_models[i];
        if(!has_class(s, "SubstitutionModel")) {
            stop("Each argument in `substitution_models` should be of class `SubstitutionModel`");
        }
    }

    List first_model = substitution_models[0];
    string model_type = get_attr(first_model, "type");

    for(ullong i = 0; i < substitution_models.size(); i++) {
        List s = substitution_models[i];
        if(get_attr(s, "type") != model_type) {
            stop("All models in `substitution_models` should have the same `type`");
        }
        if(get_attr(s, "scaling_type") != get_attr(first_model, "scaling_type")) {
            stop("All models in `substitution_models` should have the same `scaling_type`");
        }
    }
    if(!has_class(sequence, "Sequence")) {
        stop("Argument `sequence` should be of class `Sequence`");
    }
    if(get_attr(sequence, "type") != "codon") {
        stop("Argument `sequence` should be of type `codon`");
    }

    Palantir::GeneticCode g(get_genetic_code_name());

    string newick = phylogeny["newick"];
    Palantir::Phylogeny p(newick);

    mat switching_sampling = switching_model["sampling"];

    uvec states = sequence["index"];
    unsigned long long n_sites = sequence["length"];
    uvec codons = Palantir::CompoundCodon::to_codon_index(states, g);

    vector<vector<Palantir::IntervalHistory>> tree_intervals = Palantir::Simulate::switching_poisson(
        p, switching_sampling, n_sites, start_mode, switching_rate);
    List site_intervals;
    for(unsigned long long site = 0; site < n_sites; site++) {
        List intervals = interval_histories_to_list(tree_intervals[site], p);
        site_intervals.push_back(DataFrame(intervals));
        if(!compare_modes(substitution_models, intervals)) {
            stop("The modes of the `switching_model` should index the elements in `substitution_models` list");
        }
    }

    vector<vec> equilibrium;
    vector<mat> transition;
    vector<mat> sampling;

    for(ullong i = 0; i < substitution_models.size(); i++) {
        List substitution_model = substitution_models[i];
        equilibrium.push_back(substitution_model["equilibrium"]);
        transition.push_back(substitution_model["transition"]);
        sampling.push_back(substitution_model["sampling"]);
    }

    // FIX (2026-08-17): see the note in simulate_over_interval_phylogeny.
    string th_scaling_type = get_attr(List(substitution_models[0]), "scaling_type");
    vector<Palantir::SiteSimulation> sims;
    for(unsigned long long site = 0; site < n_sites; site++) {
        uvec seq(1);
        seq[0] = codons[site];

        vector<Palantir::SiteSimulation> sim = Palantir::Simulate::sequence_over_intervals(
             p, tree_intervals[site], equilibrium, transition, sampling, seq, start_mode, g, rate, segment_length, tolerance, th_scaling_type);
        if(sim.size() != 1) {
            Rcout << sim.size() << endl;
            stop("This should never happen!");
        }
        sims.push_back(sim[0]);
    }

    List substitutions = site_simulations_to_list(sims, p);

    decorate_substitutions(substitutions, model_type);

    CharacterMatrix alignment = get_alignment(sims, p, model_type);

    List simulation = List::create(
        _["phylogeny"] = phylogeny,
        _["models"] = substitution_models,
        _["substitutions"] = DataFrame::create(substitutions, _["stringsAsFactors"] = false),
        _["alignment"] = alignment,
        _["intervals"] = site_intervals,
        _["type"] = "compound_codon"
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}
