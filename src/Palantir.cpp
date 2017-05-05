#include <RcppArmadillo.h>

#include "Palantir_Core/Phylogeny.hpp"
#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/GeneticCode.hpp"
#include "Palantir_Core/SiteSimulation.hpp"
#include "Palantir_Core/Simulate.hpp"

#include "RcppPalantir.hpp"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List Phylogeny(std::string newick)
{
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

    unsigned long long size = 0;

    for(unsigned long long site = 0; site < sims.size(); site ++) {
        for(unsigned long long node = 0; node < sims[site].substitutions.size(); node ++) {
            const Palantir::SubstitutionHistory& s = sims[site].substitutions[node];
            size += s.size;
        }
    }
    arma::uvec sites(size);
    arma::uvec node_index(size);
    arma::vec time(size);
    arma::uvec state_from(size);
    arma::uvec state_to(size);

    unsigned long long i = 0;
    for(unsigned long long site = 0; site < sims.size(); site ++) {
        for(unsigned long long node = 0; node < sims[site].substitutions.size(); node ++) {
            const Palantir::SubstitutionHistory& s = sims[site].substitutions[node];
            for(unsigned long long sub = 0; sub < s.size; sub ++) {
                sites[i] = site;
                node_index[i] = node;
                time[i] = s.time[sub];
                state_from[i] = s.state_from[sub];
                state_to[i] = s.state_to[sub];
                i ++;
            }
        }
    }

    CharacterMatrix alignment = get_alignment(sims, p, substitution_model["type"]);

    List substitutions = List::create(
        _["site"] = sites,
        _["node"] = node_index,
        _["time"] = time,
        _["from"] = state_from,
        _["to"] = state_to);

    decorate_substitutions(substitutions, substitution_model["type"]);

    List simulation = List::create(
        _["phylogeny"] = tree,
        _["model"] = substitution_model,
        _["substitutions"] = DataFrame(substitutions),
        _["alignment"] = alignment,
        _["intervals"] = NULL
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}

