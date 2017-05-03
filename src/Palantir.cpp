#include <RcppArmadillo.h>

#include "Phylogeny.hpp"
#include "Palantir.hpp"
#include "HasegawaKishinoYano.hpp"
#include "MutationSelection.hpp"
#include "GeneticCode.hpp"
#include "SubstitutionHistory.hpp"
#include "Simulate.hpp"

using namespace Rcpp;
using namespace std;

std::string get_genetic_code_name() {
    Environment globals = Environment::namespace_env("PalantiR").get(".globals");
    return globals["genetic_code_name"];
}

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
        arma::vec fitness)
{
    Palantir::GeneticCode g(get_genetic_code_name());

    arma::vec equilibrium = Palantir::MutationSelection::equilibrium(population_size, mutation_rate, nucleotide_equilibrium, fitness, g);
    arma::mat transition = Palantir::MutationSelection::transition(population_size, mutation_rate, nucleotide_transition, fitness, g);
    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List ms = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["n_states"] = n_states,
        _["type"] = "codon"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
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
DataFrame simulate_over_time(
    List substitution_model,
    unsigned long long start,
    double duration,
    double rate = 1)
{
    std::string cl = substitution_model.attr("class");
    if(cl != "SubstitutionModel") {
        stop("Argument `substitution_model` should be of class `SubstitutionModel`");
    }

    arma::mat transition = substitution_model["transition"];
    arma::mat sampling = substitution_model["sampling"];

    Palantir::SubstitutionHistory sh =
        Palantir::Simulate::over_time(transition, sampling, start, duration, rate);

    Rcout << sh.size << endl;

    return DataFrame::create(
        _["time"] = sh.time,
        _["from"] = sh.state_from,
        _["to"] = sh.state_to
    );
}


// [[Rcpp::export]]
DataFrame simulate_over_phylogeny(
    List tree,
    List substitution_model,
    arma::uvec sequence,
    double rate = 1)
{
    arma::mat transition = substitution_model["transition"];
    arma::mat sampling = substitution_model["sampling"];

    string newick = tree["newick"];
    Palantir::Phylogeny p(newick);

    vector<vector<Palantir::SubstitutionHistory>> sh =
        Palantir::Simulate::sequence_over_phylogeny(
            p, transition, sampling, sequence, rate);
    Rcout << sh.size() << std::endl;
    Rcout << sh[0].size() << std::endl;

    unsigned long long size = 0;

    for(unsigned long long site = 0; site < sh.size(); site ++) {
        for(unsigned long long node = 0; node < sh[site].size(); node ++) {
            Palantir::SubstitutionHistory& s = sh[site][node];
            size += s.size;
        }
    }
    Rcout << size << std::endl;
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
                time[i] = s.time[i];
                state_from[i] = s.state_from[i];
                state_to[i] = s.state_to[i];
                i ++;
            }
        }
    }

    return DataFrame::create(
        _["site"] = sites,
        _["node"] = node_index,
        _["time"] = time,
        _["from"] = state_from,
        _["to"] = state_to);
}
