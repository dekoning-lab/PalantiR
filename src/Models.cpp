#include <RcppArmadillo.h>

#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/HasegawaKishinoYano.hpp"
#include "Palantir_Core/GeneralTimeReversible.hpp"
#include "Palantir_Core/MutationSelection.hpp"
#include "Palantir_Core/CoEvolution.hpp"
#include "Palantir_Core/MarkovModel.hpp"
#include "Palantir_Core/MarkovModulated.hpp"

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
        _["transition_rate"] = transition_rate,
        _["transversion_rate"] = transversion_rate,
        _["n_states"] = n_states,
        _["type"] = "nucleotide"
    );

    hky.attr("class") = "SubstitutionModel";
    return hky;
}

//[[Rcpp::export]]
List GeneralTimeReversible(arma::vec equilibrium, arma::mat exchangeability)
{
    arma::mat transition = Palantir::GeneralTimeReversible::transition(equilibrium, exchangeability);
    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List gtr = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["exchangeability"] = exchangeability,
        _["n_states"] = n_states,
        _["type"] = "exchangeable"
    );

    gtr.attr("class") = "SubstitutionModel";
    return gtr;
}

// [[Rcpp::export]]
List MutationSelection(
        unsigned long long population_size,
        double mutation_rate,
        List nucleotide_model,
        arma::vec fitness,
        std::string scaling_type = "synonymous")
{
    Palantir::GeneticCode g(get_genetic_code_name());

    if(!has_class(nucleotide_model, "SubstitutionModel") || get_attr(nucleotide_model, "type") != "nucleotide") {
        stop("Argument `nucleotide_model` should be a nucleotide substitution model");
    }

    arma::vec nucleotide_equilibrium = nucleotide_model["equilibrium"];
    arma::mat nucleotide_transition = nucleotide_model["transition"];

    arma::vec equilibrium = Palantir::MutationSelection::equilibrium(
        population_size, mutation_rate, nucleotide_equilibrium, fitness, g);
    arma::mat transition = Palantir::MutationSelection::transition(
        population_size, mutation_rate, nucleotide_transition, fitness, g);
    double scaling = Palantir::MutationSelection::scaling(
        equilibrium, transition, scaling_type, g);
    transition /= scaling;

    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List ms = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["population_size"] = population_size,
        _["fitness"] = fitness,
        _["nucleotide_model"] = nucleotide_model,
        _["scaling"] = scaling,
        _["scaling_type"] = scaling_type,
        _["n_states"] = n_states,
        _["type"] = "codon"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}

//[[Rcpp::export]]
List CoEvolution(
    unsigned long long population_size,
    double mutation_rate,
    List nucleotide_model,
    arma::vec fitness_1,
    arma::vec fitness_2,
    arma::mat delta,
    std::string scaling_type = "synonymous")
{
    Palantir::GeneticCode g(get_genetic_code_name());

    if(!has_class(nucleotide_model, "SubstitutionModel") || get_attr(nucleotide_model, "type") != "nucleotide") {
        stop("Argument `nucleotide_model` should be a nucleotide substitution model");
    }
    arma::vec nucleotide_equilibrium = nucleotide_model["equilibrium"];
    arma::mat nucleotide_transition = nucleotide_model["transition"];

    arma::vec equilibrium = Palantir::CoEvolution::equilibrium(
        population_size, mutation_rate, nucleotide_equilibrium, fitness_1, fitness_2, delta, g);
    arma::mat transition = Palantir::CoEvolution::transition(
        population_size, mutation_rate, nucleotide_transition, fitness_1, fitness_2, delta, g);
    double scaling = Palantir::CoEvolution::scaling(equilibrium, transition, scaling_type, g);
    transition /= scaling;

    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List ms = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["population_size"] = population_size,
        _["fitness_1"] = fitness_1,
        _["fitness_2"] = fitness_2,
        _["nucleotide_model"] = nucleotide_model,
        _["scaling"] = scaling,
        _["scaling_type"] = scaling_type,
        _["n_states"] = n_states,
        _["type"] = "codon_pair"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}

//[[Rcpp::export]]
List MarkovModulatedMutationSelection(
    List mutation_selection_models,
    List switching_model)
{
    Palantir::GeneticCode g(get_genetic_code_name());
    if(!has_class(switching_model, "SubstitutionModel") || get_attr(switching_model, "type") != "exchangeable") {
        stop("Argument `switching_model` should be an exchangeable substitution model");
    }
    vec switching_equilibrium = switching_model["equilibrium"];
    mat exchangeability = switching_model["exchangeability"];
    vector<const vec> substitution_equilibrium;
    vector<const mat> substitution_transition;

    for(ullong i = 0; i < mutation_selection_models.size(); i++) {
        List ms_model = mutation_selection_models[i];
        substitution_equilibrium.push_back(ms_model["equilibrium"]);
        substitution_transition.push_back(ms_model["transition"]);
    }

    mat transition = Palantir::MarkovModulated::transition(
        substitution_transition, exchangeability, switching_equilibrium, g);

    vec equilibrium = Palantir::MarkovModel::solve_equilibrium(transition);

    unsigned long long n_states = equilibrium.n_elem;

    mat sampling = Palantir::sampling(transition);

    List ms = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["mutation_selection_models"] = mutation_selection_models,
        _["n_states"] = n_states,
        _["type"] = "compound_codon"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}
