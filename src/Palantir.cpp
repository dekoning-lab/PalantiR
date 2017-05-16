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
    List tree,
    List substitution_model,
    List sequence,
    double rate = 1)
{
    // Type checking
    if(!has_class(tree, "Phylogeny")) {
        stop("Argument `tree` should be of class `Phylogeny`");
    }
    if(!has_class(substitution_model, "SubstitutionModel")) {
        stop("Argument `substitution_model` should be of class `SubstitutionModel`");
    }
    if(!has_class(sequence, "Sequence")) {
        stop("Argument `sequence` should be of class `Sequence`");
    }

    arma::mat transition = substitution_model["transition"];
    arma::mat sampling = substitution_model["sampling"];

    string newick = tree["newick"];
    Palantir::Phylogeny p(newick);
    uvec states = sequence["index"];

    vector<Palantir::SiteSimulation> sims =
        Palantir::Simulate::sequence_over_phylogeny(
            p, transition, sampling, states, rate);

    List substitutions = site_simulations_to_list(sims, p);
    decorate_substitutions(substitutions, substitution_model["type"]);

    CharacterMatrix alignment = get_alignment(sims, p, substitution_model["type"]);

    List simulation = List::create(
        _["phylogeny"] = tree,
        _["model"] = substitution_model,
        _["substitutions"] = DataFrame(substitutions),
        _["alignment"] = alignment,
        _["intervals"] = NULL,
        _["type"] = substitution_model["type"]
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}

//[[Rcpp::export]]
List simulate_over_interval_phylogeny(
    List tree,
    List mode_tree,
    List substitution_models,
    List sequence,
    double rate = 1,
    double segment_length = 0.001,
    double tolerance = 0.001)
{

    //Type checking
    List first_model = substitution_models[0];
    string model_type = get_attr(first_model, "type");

    if(!has_class(tree, "Phylogeny")) {
        stop("Argument `tree` should be of class `Phylogeny`");
    }
    if(!has_class(mode_tree, "Phylogeny")) {
        stop("Argument `mode_tree` should be of class `Phylogeny`");
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
    if(get_attr(sequence, "type") != "compound_codon") {
        stop("Argument `sequence` should be of type `compound_codon`");
    }

    Palantir::GeneticCode g(get_genetic_code_name());

    string newick = tree["newick"];
    Palantir::Phylogeny p(newick);

    string mode_newick = mode_tree["newick"];
    Palantir::Phylogeny pa(mode_newick);

    vector<Palantir::IntervalHistory> tree_intervals = p.to_intervals(pa);
    List intervals = interval_histories_to_list(tree_intervals, p);

    if(!compare_modes(substitution_models, intervals)) {
        stop("The modes on the `mode_tree` phylogeny should index the elements in `substitution_models` list");
    }

    uvec states = sequence["index"];
    uvec modes = Palantir::CompoundCodon::to_mode_index(states, g);
    uvec codons = Palantir::CompoundCodon::to_codon_index(states, g);
    ullong start_mode = modes[0];

    if(any(modes != start_mode)) {
        stop("Every site in the sequence should have the same mode");
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
        _["phylogeny"] = tree,
        _["models"] = substitution_models,
        _["substitutions"] = DataFrame(substitutions),
        _["alignment"] = alignment,
        _["intervals"] = DataFrame(intervals),
        _["type"] = "compound_codon"
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}

//[[Rcpp::export]]
List simulate_with_nested_heterogeneity(
    List tree,
    List switching_model,
    List substitution_models,
    List sequence,
    double rate = 1,
    double switching_rate = 1,
    double segment_length = 0.001,
    double tolerance = 0.001)
{
    //Type checking
    List first_model = substitution_models[0];
    string model_type = get_attr(first_model, "type");

    if(!has_class(tree, "Phylogeny")) {
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
    if(get_attr(sequence, "type") != "compound_codon") {
        stop("Argument `sequence` should be of type `compound_codon`");
    }

    Palantir::GeneticCode g(get_genetic_code_name());

    string newick = tree["newick"];
    Palantir::Phylogeny p(newick);

    mat switching_transition = switching_model["transition"];
    mat switching_sampling = switching_model["sampling"];

    uvec states = sequence["index"];
    uvec modes = Palantir::CompoundCodon::to_mode_index(states, g);
    uvec codons = Palantir::CompoundCodon::to_codon_index(states, g);
    ullong start_mode = modes[0];
    if(any(modes != start_mode)) {
        stop("Every site in the sequence should have the same mode");
    }

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
        _["phylogeny"] = tree,
        _["models"] = substitution_models,
        _["substitutions"] = DataFrame(substitutions),
        _["alignment"] = alignment,
        _["intervals"] = DataFrame(intervals),
        _["type"] = "compound_codon"
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}
