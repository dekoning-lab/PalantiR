#include <RcppArmadillo.h>

#include "Phylogeny.hpp"
#include "Palantir.hpp"
#include "HasegawaKishinoYano.hpp"
#include "MutationSelection.hpp"
#include "GeneticCode.hpp"
#include "SubstitutionHistory.hpp"
#include "Simulate.hpp"

#include "RcppUtil.hpp"

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
List HasegawaKishinoYano(arma::vec equilibrium, double transition_rate = 1, double transversion_rate = 1)
{
    arma::mat transition = Palantir::HasegawaKishinoYano::transition(equilibrium, transition_rate, transversion_rate);
    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List hky = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["n_states"] = n_states,
        _["type"] = "nucleotide"
    );

    hky.attr("class") = "SubstitutionModel";
    return hky;
}

// [[Rcpp::export]]
List MutationSelection(
        unsigned long long population_size,
        double mutation_rate,
        arma::vec nucleotide_equilibrium,
        arma::mat nucleotide_transition,
        arma::vec fitness,
        std::string scaling_type = "synonymous")
{
    Palantir::GeneticCode g(get_genetic_code_name());

    arma::vec equilibrium = Palantir::MutationSelection::equilibrium(population_size, mutation_rate, nucleotide_equilibrium, fitness, g);
    arma::mat transition = Palantir::MutationSelection::transition(population_size, mutation_rate, nucleotide_transition, fitness, g);
    double scaling = Palantir::MutationSelection::scaling(equilibrium, transition, scaling_type, g);
    transition /= scaling;

    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List ms = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["scaling"] = scaling,
        _["n_states"] = n_states,
        _["type"] = "codon"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}

CharacterVector predicate(uvec predicate, string yes, string no)
{
    CharacterVector c(predicate.size());
    for(unsigned long long i = 0; i < predicate.size(); i++) {
        if(predicate[i]) {
            c[i] = yes;
        } else {
            c[i] = no;
        }
    }
    return c;
}

// [[Rcpp::export]]
DataFrame decorate_codon_substitutions(DataFrame substitutions)
{
    Palantir::GeneticCode g(get_genetic_code_name());
    arma::uvec state_from = substitutions["from"];
    arma::uvec state_to = substitutions["to"];

    substitutions["codon_from"] = Palantir::Codon::to_string(state_from, g);
    substitutions["codon_to"] = Palantir::Codon::to_string(state_to, g);
    substitutions["amino_acid_from"] = Palantir::Codon::to_amino_acid(state_from, g);
    substitutions["amino_acid_to"] = Palantir::Codon::to_amino_acid(state_to, g);

    arma::uvec synonymous = Palantir::Codon::synonymous(state_from, state_to, g);
    substitutions["synonymous"] = LogicalVector(synonymous.begin(), synonymous.end());
    substitutions["color"] = predicate(synonymous, "#16A35E", "#E84B2A");

    return substitutions;
}

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
    arma::uvec sequence,
    double rate = 1)
{
    if(!has_class(tree, "Phylogeny")) {
        stop("Argument `tree` should be of class `Phylogeny`");
    }
    if(!has_class(substitution_model, "SubstitutionModel")) {
        stop("Argument `substitution_model` should be of class `SubstitutionModel`");
    }

    arma::mat transition = substitution_model["transition"];
    arma::mat sampling = substitution_model["sampling"];

    string newick = tree["newick"];
    Palantir::Phylogeny p(newick);

    vector<vector<Palantir::SubstitutionHistory>> sh =
        Palantir::Simulate::sequence_over_phylogeny(
            p, transition, sampling, sequence, rate);

    unsigned long long size = 0;

    for(unsigned long long site = 0; site < sh.size(); site ++) {
        for(unsigned long long node = 0; node < sh[site].size(); node ++) {
            Palantir::SubstitutionHistory& s = sh[site][node];
            size += s.size;
        }
    }
    arma::uvec sites(size);
    arma::uvec node_index(size);
    arma::vec time(size);
    arma::uvec state_from(size);
    arma::uvec state_to(size);

    unsigned long long i = 0;
    for(unsigned long long site = 0; site < sh.size(); site ++) {
        for(unsigned long long node = 0; node < sh[site].size(); node ++) {
            Palantir::SubstitutionHistory& s = sh[site][node];
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

    DataFrame substitutions = DataFrame::create(
        _["site"] = sites,
        _["node"] = node_index,
        _["time"] = time,
        _["from"] = state_from,
        _["to"] = state_to);


    List simulation = List::create(
        _["phylogeny"] = tree,
        _["model"] = substitution_model,
        _["substitutions"] = substitutions,
        _["intervals"] = NULL
    );

    simulation.attr("class") = "Simulation";

    return simulation;
}
