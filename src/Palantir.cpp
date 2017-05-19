#include <RcppArmadillo.h>

#include "Palantir_Core/Phylogeny.hpp"
#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/GeneticCode.hpp"
#include "Palantir_Core/SiteSimulation.hpp"
#include "Palantir_Core/Simulate.hpp"
#include "Palantir_Core/Util.hpp"

#include "RcppPalantir.hpp"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List Phylogeny(std::string newick_path)
{
    string newick;

    ifstream f(newick_path);
    if(f.good()) {
        newick = file_to_string(newick_path);
    } else {
        stop("Could not find file " + newick_path);
    }
    //TODO: need error checking and file parsing
    Palantir::Phylogeny tree(newick);

    List phylo = List::create(
        _["json"] = tree.to_JSON(),
        _["string"] = tree.to_string(),
        _["newick"] = newick,
        _["n_nodes"] = tree.n_nodes
    );
    phylo.attr("class") = "Phylogeny";
    return phylo;
}

// [[Rcpp::export]]
arma::vec equilibrium_to_fitness(arma::vec equilibrium, unsigned long long population_size)
{
    return Palantir::equilibrium_to_fitness(equilibrium, population_size);
}

// [[Rcpp::export]]
List simulate_over_phylogeny(
    List phylogeny,
    List model,
    List sequence,
    double rate = 1)
{
    // Type checking
    if(!has_class(phylogeny, "Phylogeny")) {
        stop("Argument `tree` should be of class `Phylogeny`");
    }
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
        _["intervals"] = NULL,
        _["type"] = model["type"]
    );

    simulation.attr("class") = "Simulation";

    return simulation;
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
    double tolerance = 0.001)
{

    //Type checking
    List first_model = models[0];
    string model_type = get_attr(first_model, "type");

    if(!has_class(phylogeny, "Phylogeny")) {
        stop("Argument `tree` should be of class `Phylogeny`");
    }
    if(!has_class(mode_phylogeny, "Phylogeny")) {
        stop("Argument `mode_tree` should be of class `Phylogeny`");
    }
    for(ullong i = 0; i < models.size(); i++) {
        List s = models[i];
        if(!has_class(s, "SubstitutionModel")) {
            stop("Each argument in `substitution_models` should be of class `SubstitutionModel`");
        }

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
    if(get_attr(sequence, "type") != "compound_codon") {
        stop("Argument `sequence` should be of type `compound_codon`");
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

    vector<Palantir::SiteSimulation> sims = Palantir::Simulate::sequence_over_intervals(
        p, tree_intervals, equilibrium, transition, sampling, codons, start_mode, g, rate, segment_length, tolerance);

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
List simulate_with_nested_heterogeneity(
    List phylogeny,
    List switching_model,
    List substitution_models,
    List sequence,
    unsigned long long start_mode,
    double rate = 1,
    double switching_rate = 1,
    double segment_length = 0.001,
    double tolerance = 0.001)
{
    //Type checking
    List first_model = substitution_models[0];
    string model_type = get_attr(first_model, "type");

    if(!has_class(phylogeny, "Phylogeny")) {
        stop("Argument `tree` should be of class `Phylogeny`");
    }
    for(ullong i = 0; i < substitution_models.size(); i++) {
        List s = substitution_models[i];
        if(!has_class(s, "SubstitutionModel")) {
            stop("Each argument in `substitution_models` should be of class `SubstitutionModel`");
        }

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

    vector<Palantir::SiteSimulation> sims = Palantir::Simulate::sequence_over_intervals(
        p, tree_intervals, equilibrium, transition, sampling, codons, start_mode, g, rate, segment_length, tolerance);

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
List simulate_with_poisson_heterogeneity(
    List phylogeny,
    List switching_model,
    List substitution_models,
    List sequence,
    unsigned long long start_mode,
    double rate = 1,
    double switching_rate = 1,
    double segment_length = 0.001,
    double tolerance = 0.001)
{
    //Type checking
    List first_model = substitution_models[0];
    string model_type = get_attr(first_model, "type");

    if(!has_class(phylogeny, "Phylogeny")) {
        stop("Argument `tree` should be of class `Phylogeny`");
    }
    for(ullong i = 0; i < substitution_models.size(); i++) {
        List s = substitution_models[i];
        if(!has_class(s, "SubstitutionModel")) {
            stop("Each argument in `substitution_models` should be of class `SubstitutionModel`");
        }

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

    vector<Palantir::SiteSimulation> sims;
    for(unsigned long long site = 0; site < n_sites; site++) {
        uvec seq(1);
        seq[0] = codons[site];

        vector<Palantir::SiteSimulation> sim = Palantir::Simulate::sequence_over_intervals(
             p, tree_intervals[site], equilibrium, transition, sampling, seq, start_mode, g, rate, segment_length, tolerance);
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
