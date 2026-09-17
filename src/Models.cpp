#include <RcppArmadillo.h>

#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/HasegawaKishinoYano.hpp"
#include "Palantir_Core/GeneralTimeReversible.hpp"
#include "Palantir_Core/MutationSelection.hpp"
#include "Palantir_Core/CodonModel.hpp"
#include "Palantir_Core/GoldmanYang94.hpp"
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
        double population_size,
        double mutation_rate,
        List nucleotide_model,
        arma::vec fitness,
        std::string scaling_type = "synonymous")
{
    Palantir::GeneticCode g(get_genetic_code_name());
    scaling_type = Palantir::CodonModel::canonical_scaling_type(scaling_type);

    if(!has_class(nucleotide_model, "SubstitutionModel") || get_attr(nucleotide_model, "type") != "nucleotide") {
        stop("Argument `nucleotide_model` should be a nucleotide substitution model");
    }

    // FIX (2026-08-20, M3/C8): validate before anything is constructed. Ne = 0
    // used to produce an all-NaN model that printed normally, a negative Ne
    // wrapped to a huge unsigned value, and a mis-sized `fitness` was read out
    // of bounds by armadillo's unchecked accessors.
    unsigned long long N = as_count(population_size, "population_size", 1);
    as_positive_rate(mutation_rate, "mutation_rate");
    check_length(fitness.n_elem, Palantir::AminoAcid::size, "fitness");

    arma::vec nucleotide_equilibrium = nucleotide_model["equilibrium"];
    arma::mat nucleotide_transition = nucleotide_model["transition"];

    arma::vec equilibrium = Palantir::MutationSelection::equilibrium(
        N, mutation_rate, nucleotide_equilibrium, fitness, g);
    arma::mat transition = Palantir::MutationSelection::transition(
        N, mutation_rate, nucleotide_transition, fitness, g);
    double neutral_total_rate = NA_REAL;
    double neutral_synonymous_rate = NA_REAL;
    double synonymous_opportunities = NA_REAL;
    double scaling;
    if(scaling_type == "dS") {
        // A dS branch is neutral nucleotide time, not one realised
        // synonymous event per selected codon. Build the mutation-only
        // reference once for this constructor. It is independent of the
        // supplied amino-acid fitness profile, so site-heterogeneous models
        // that share a mutation process retain a common clock.
        arma::vec neutral_fitness(Palantir::AminoAcid::size, fill::zeros);
        arma::vec neutral_equilibrium = Palantir::MutationSelection::equilibrium(
            N, mutation_rate, nucleotide_equilibrium, neutral_fitness, g);
        arma::mat neutral_transition = Palantir::MutationSelection::transition(
            N, mutation_rate, nucleotide_transition, neutral_fitness, g);
        scaling = Palantir::CodonModel::neutral_dS_scaling(
            neutral_equilibrium, neutral_transition, g);
        neutral_total_rate = 3.0 * scaling;
        neutral_synonymous_rate = Palantir::CodonModel::scaling(
            neutral_equilibrium, neutral_transition,
            "synonymous-per-codon", g);
        synonymous_opportunities =
            3.0 * neutral_synonymous_rate / neutral_total_rate;
    } else {
        scaling = Palantir::MutationSelection::scaling(
            equilibrium, transition, scaling_type, g);
    }
    transition /= scaling;

    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List ms = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["population_size"] = N,
        _["fitness"] = fitness,
        _["nucleotide_model"] = nucleotide_model,
        _["scaling"] = scaling,
        _["scaling_type"] = scaling_type,
        _["neutral_total_rate"] = neutral_total_rate,
        _["neutral_synonymous_rate"] = neutral_synonymous_rate,
        _["synonymous_opportunities"] = synonymous_opportunities,
        _["n_states"] = n_states,
        _["type"] = "codon"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}

//[[Rcpp::export]]
List CoEvolution(
    double population_size,
    double mutation_rate,
    List nucleotide_model,
    arma::vec fitness_1,
    arma::vec fitness_2,
    arma::mat delta,
    std::string scaling_type = "synonymous")
{
    Palantir::GeneticCode g(get_genetic_code_name());
    scaling_type = Palantir::CodonModel::canonical_scaling_type(scaling_type);

    if(!has_class(nucleotide_model, "SubstitutionModel") || get_attr(nucleotide_model, "type") != "nucleotide") {
        stop("Argument `nucleotide_model` should be a nucleotide substitution model");
    }

    // FIX (2026-08-20, M3/C8): see the note in MutationSelection above. The
    // coevolution kernel indexes fitness_1/fitness_2 by amino acid and delta by
    // the amino-acid pair, all with unchecked accessors: a 2x2 delta used to
    // give a "model" whose equilibrium summed to 24.
    unsigned long long N = as_count(population_size, "population_size", 1);
    as_positive_rate(mutation_rate, "mutation_rate");
    check_length(fitness_1.n_elem, Palantir::AminoAcid::size, "fitness_1");
    check_length(fitness_2.n_elem, Palantir::AminoAcid::size, "fitness_2");
    check_square(delta.n_rows, delta.n_cols, Palantir::AminoAcid::size, "delta");

    arma::vec nucleotide_equilibrium = nucleotide_model["equilibrium"];
    arma::mat nucleotide_transition = nucleotide_model["transition"];

    arma::vec equilibrium = Palantir::CoEvolution::equilibrium(
        N, mutation_rate, nucleotide_equilibrium, fitness_1, fitness_2, delta, g);
    arma::mat transition = Palantir::CoEvolution::transition(
        N, mutation_rate, nucleotide_transition, fitness_1, fitness_2, delta, g);
    double neutral_total_rate = NA_REAL;
    double neutral_synonymous_rate = NA_REAL;
    double synonymous_opportunities = NA_REAL;
    double scaling;
    if(scaling_type == "dS") {
        // The neutral reference for a codon pair is the Kronecker sum of two
        // identical single-codon mutation processes. Its rate per codon is
        // therefore exactly the single-codon rate. Compute that 61-state
        // reference directly instead of allocating a second 3721 x 3721
        // matrix; this keeps dS construction negligible beside CoEvolution's
        // selected generator.
        arma::vec neutral_fitness(Palantir::AminoAcid::size, fill::zeros);
        arma::vec neutral_equilibrium = Palantir::MutationSelection::equilibrium(
            N, mutation_rate, nucleotide_equilibrium, neutral_fitness, g);
        arma::mat neutral_transition = Palantir::MutationSelection::transition(
            N, mutation_rate, nucleotide_transition, neutral_fitness, g);
        scaling = Palantir::CodonModel::neutral_dS_scaling(
            neutral_equilibrium, neutral_transition, g);
        neutral_total_rate = 3.0 * scaling;
        neutral_synonymous_rate = Palantir::CodonModel::scaling(
            neutral_equilibrium, neutral_transition,
            "synonymous-per-codon", g);
        synonymous_opportunities =
            3.0 * neutral_synonymous_rate / neutral_total_rate;
    } else {
        scaling = Palantir::CoEvolution::scaling(
            equilibrium, transition, scaling_type, g);
    }
    transition /= scaling;

    arma::mat sampling = Palantir::sampling(transition);
    unsigned long long n_states = equilibrium.n_elem;

    List ms = List::create(
        _["equilibrium"] = equilibrium,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["population_size"] = N,
        _["fitness_1"] = fitness_1,
        _["fitness_2"] = fitness_2,
        _["nucleotide_model"] = nucleotide_model,
        _["scaling"] = scaling,
        _["scaling_type"] = scaling_type,
        _["neutral_total_rate"] = neutral_total_rate,
        _["neutral_synonymous_rate"] = neutral_synonymous_rate,
        _["synonymous_opportunities"] = synonymous_opportunities,
        _["n_states"] = n_states,
        _["type"] = "codon_pair"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}

// [[Rcpp::export(name = ".GoldmanYang94Cpp")]]
List GoldmanYang94Cpp(
        arma::vec equilibrium,
        double omega,
        double kappa,
        std::string frequency_model,
        std::string scaling_type = "substitution")
{
    Palantir::GeneticCode g(get_genetic_code_name());
    scaling_type = Palantir::CodonModel::canonical_scaling_type(scaling_type);

    if(!std::isfinite(omega) || omega <= 0) {
        stop("Argument `omega` should be a finite number greater than 0");
    }
    if(!std::isfinite(kappa) || kappa <= 0) {
        stop("Argument `kappa` should be a finite number greater than 0");
    }
    if(equilibrium.n_elem != g.size) {
        stop("Argument `equilibrium` should have " +
             std::to_string(g.size) + " elements (one per sense codon under "
             "the active genetic code), not " +
             std::to_string(equilibrium.n_elem));
    }
    if(!equilibrium.is_finite() || equilibrium.min() <= 0) {
        stop("Argument `equilibrium` should contain strictly positive finite "
             "sense-codon frequencies");
    }
    double equilibrium_sum = arma::sum(equilibrium);
    if(std::abs(equilibrium_sum - 1.0) > 1e-10) {
        stop("Argument `equilibrium` should sum to 1 (received " +
             std::to_string(equilibrium_sum) + ")");
    }

    if(frequency_model == "FEqual") {
        frequency_model = "Fequal";
    }
    if(frequency_model != "Fequal" && frequency_model != "F1x4" &&
       frequency_model != "F3x4" && frequency_model != "F61") {
        stop("Argument `frequency_model` should be one of \"Fequal\", "
             "\"F1x4\", \"F3x4\" or \"F61\"");
    }

    arma::mat transition = Palantir::GoldmanYang94::transition(
        equilibrium, omega, kappa, g);
    double neutral_total_rate = NA_REAL;
    double neutral_synonymous_rate = NA_REAL;
    double synonymous_opportunities = NA_REAL;
    double scaling;
    if(scaling_type == "dS") {
        // Define the clock from omega=1 so changing selection does not change
        // the meaning of the supplied branch length.
        arma::mat neutral_transition = Palantir::GoldmanYang94::transition(
            equilibrium, 1.0, kappa, g);
        scaling = Palantir::CodonModel::neutral_dS_scaling(
            equilibrium, neutral_transition, g);
        neutral_total_rate = 3.0 * scaling;
        neutral_synonymous_rate = Palantir::CodonModel::scaling(
            equilibrium, neutral_transition, "synonymous-per-codon", g);
        synonymous_opportunities =
            3.0 * neutral_synonymous_rate / neutral_total_rate;
    } else {
        scaling = Palantir::CodonModel::scaling(
            equilibrium, transition, scaling_type, g);
    }
    if(!std::isfinite(scaling) || scaling <= 0) {
        stop("The requested scaling class has a non-positive or non-finite "
             "stationary rate for this GY94 model");
    }
    transition /= scaling;
    arma::mat sampling = Palantir::sampling(transition);

    CharacterVector codon_names(g.size);
    for(ullong i = 0; i < g.size; i++) {
        codon_names[i] = g.sequence[i];
    }
    NumericVector codon_frequencies = wrap(equilibrium);
    codon_frequencies.attr("names") = codon_names;

    List gy = List::create(
        _["model"] = "GY94",
        _["genetic_code"] = get_genetic_code_name(),
        _["equilibrium"] = codon_frequencies,
        _["codon_frequencies"] = codon_frequencies,
        _["transition"] = transition,
        _["sampling"] = sampling,
        _["omega"] = omega,
        _["kappa"] = kappa,
        _["frequency_model"] = frequency_model,
        _["scaling"] = scaling,
        _["scaling_type"] = scaling_type,
        _["neutral_total_rate"] = neutral_total_rate,
        _["neutral_synonymous_rate"] = neutral_synonymous_rate,
        _["synonymous_opportunities"] = synonymous_opportunities,
        _["n_states"] = g.size,
        _["type"] = "codon"
    );
    gy.attr("class") = CharacterVector::create(
        "GoldmanYang94", "SubstitutionModel");
    return gy;
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
    // FIX (2026-08-20, MIN7): an empty or mistyped model list used to fail deep
    // inside the kernel (or not at all) instead of at the boundary.
    if(mutation_selection_models.size() == 0) {
        stop("Argument `mutation_selection_models` should contain at least one model");
    }

    vec switching_equilibrium = switching_model["equilibrium"];
    mat exchangeability = switching_model["exchangeability"];
    vector<vec> substitution_equilibrium;
    vector<mat> substitution_transition;
    string scaling_type;

    for(ullong i = 0; i < mutation_selection_models.size(); i++) {
        List ms_model = mutation_selection_models[i];
        if(!has_class(ms_model, "SubstitutionModel")) {
            stop("Each argument in `mutation_selection_models` should be of class `SubstitutionModel`");
        }
        if(get_attr(ms_model, "type") != "codon") {
            stop("Each argument in `mutation_selection_models` should be a single-codon model");
        }
        string component_scaling = get_attr(ms_model, "scaling_type");
        if(i == 0) {
            scaling_type = component_scaling;
        } else if(component_scaling != scaling_type) {
            stop("All `mutation_selection_models` should use the same `scaling_type`");
        }
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
        _["scaling_type"] = scaling_type,
        _["n_states"] = n_states,
        _["type"] = "compound_codon"
    );

    ms.attr("class") = "SubstitutionModel";
    return ms;
}
