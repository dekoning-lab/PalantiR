#include <RcppArmadillo.h>

#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/HasegawaKishinoYano.hpp"
#include "Palantir_Core/MutationSelection.hpp"

#include "RcppPalantir.hpp"

using namespace Rcpp;

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
