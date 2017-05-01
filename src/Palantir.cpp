#include <RcppArmadillo.h>

#include "Phylogeny.hpp"
#include "Palantir.hpp"
#include "HasegawaKishinoYano.hpp"
#include "MutationSelection.hpp"
#include "GeneticCode.hpp"

using namespace Rcpp;

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
        _["n_states"] = n_states
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
        _["n_states"] = n_states
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}
