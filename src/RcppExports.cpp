// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GeneticCode
StringVector GeneticCode();
RcppExport SEXP _PalantiR_GeneticCode() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(GeneticCode());
    return rcpp_result_gen;
END_RCPP
}
// HasegawaKishinoYano
List HasegawaKishinoYano(arma::vec equilibrium, double transition_rate, double transversion_rate);
RcppExport SEXP _PalantiR_HasegawaKishinoYano(SEXP equilibriumSEXP, SEXP transition_rateSEXP, SEXP transversion_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type equilibrium(equilibriumSEXP);
    Rcpp::traits::input_parameter< double >::type transition_rate(transition_rateSEXP);
    Rcpp::traits::input_parameter< double >::type transversion_rate(transversion_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(HasegawaKishinoYano(equilibrium, transition_rate, transversion_rate));
    return rcpp_result_gen;
END_RCPP
}
// GeneralTimeReversible
List GeneralTimeReversible(arma::vec equilibrium, arma::mat exchangeability);
RcppExport SEXP _PalantiR_GeneralTimeReversible(SEXP equilibriumSEXP, SEXP exchangeabilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type equilibrium(equilibriumSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type exchangeability(exchangeabilitySEXP);
    rcpp_result_gen = Rcpp::wrap(GeneralTimeReversible(equilibrium, exchangeability));
    return rcpp_result_gen;
END_RCPP
}
// MutationSelection
List MutationSelection(unsigned long long population_size, double mutation_rate, List nucleotide_model, arma::vec fitness, std::string scaling_type);
RcppExport SEXP _PalantiR_MutationSelection(SEXP population_sizeSEXP, SEXP mutation_rateSEXP, SEXP nucleotide_modelSEXP, SEXP fitnessSEXP, SEXP scaling_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned long long >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< List >::type nucleotide_model(nucleotide_modelSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fitness(fitnessSEXP);
    Rcpp::traits::input_parameter< std::string >::type scaling_type(scaling_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(MutationSelection(population_size, mutation_rate, nucleotide_model, fitness, scaling_type));
    return rcpp_result_gen;
END_RCPP
}
// CoEvolution
List CoEvolution(unsigned long long population_size, double mutation_rate, List nucleotide_model, arma::vec fitness_1, arma::vec fitness_2, arma::mat delta, std::string scaling_type);
RcppExport SEXP _PalantiR_CoEvolution(SEXP population_sizeSEXP, SEXP mutation_rateSEXP, SEXP nucleotide_modelSEXP, SEXP fitness_1SEXP, SEXP fitness_2SEXP, SEXP deltaSEXP, SEXP scaling_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned long long >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< List >::type nucleotide_model(nucleotide_modelSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fitness_1(fitness_1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fitness_2(fitness_2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< std::string >::type scaling_type(scaling_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(CoEvolution(population_size, mutation_rate, nucleotide_model, fitness_1, fitness_2, delta, scaling_type));
    return rcpp_result_gen;
END_RCPP
}
// MarkovModulatedMutationSelection
List MarkovModulatedMutationSelection(List mutation_selection_models, List switching_model);
RcppExport SEXP _PalantiR_MarkovModulatedMutationSelection(SEXP mutation_selection_modelsSEXP, SEXP switching_modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type mutation_selection_models(mutation_selection_modelsSEXP);
    Rcpp::traits::input_parameter< List >::type switching_model(switching_modelSEXP);
    rcpp_result_gen = Rcpp::wrap(MarkovModulatedMutationSelection(mutation_selection_models, switching_model));
    return rcpp_result_gen;
END_RCPP
}
// Phylogeny
List Phylogeny(std::string newick_path, std::string type);
RcppExport SEXP _PalantiR_Phylogeny(SEXP newick_pathSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type newick_path(newick_pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Phylogeny(newick_path, type));
    return rcpp_result_gen;
END_RCPP
}
// equilibrium_to_fitness
arma::vec equilibrium_to_fitness(arma::vec equilibrium, unsigned long long population_size);
RcppExport SEXP _PalantiR_equilibrium_to_fitness(SEXP equilibriumSEXP, SEXP population_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type equilibrium(equilibriumSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type population_size(population_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(equilibrium_to_fitness(equilibrium, population_size));
    return rcpp_result_gen;
END_RCPP
}
// simulate_over_phylogeny
List simulate_over_phylogeny(List phylogeny, List model, List sequence, double rate);
RcppExport SEXP _PalantiR_simulate_over_phylogeny(SEXP phylogenySEXP, SEXP modelSEXP, SEXP sequenceSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type phylogeny(phylogenySEXP);
    Rcpp::traits::input_parameter< List >::type model(modelSEXP);
    Rcpp::traits::input_parameter< List >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_over_phylogeny(phylogeny, model, sequence, rate));
    return rcpp_result_gen;
END_RCPP
}
// phylogeny_to_intervals
DataFrame phylogeny_to_intervals(List phylogeny, List mode_phylogeny);
RcppExport SEXP _PalantiR_phylogeny_to_intervals(SEXP phylogenySEXP, SEXP mode_phylogenySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type phylogeny(phylogenySEXP);
    Rcpp::traits::input_parameter< List >::type mode_phylogeny(mode_phylogenySEXP);
    rcpp_result_gen = Rcpp::wrap(phylogeny_to_intervals(phylogeny, mode_phylogeny));
    return rcpp_result_gen;
END_RCPP
}
// simulate_over_interval_phylogeny
List simulate_over_interval_phylogeny(List phylogeny, List mode_phylogeny, List models, List sequence, unsigned long long start_mode, double rate, double segment_length, double tolerance);
RcppExport SEXP _PalantiR_simulate_over_interval_phylogeny(SEXP phylogenySEXP, SEXP mode_phylogenySEXP, SEXP modelsSEXP, SEXP sequenceSEXP, SEXP start_modeSEXP, SEXP rateSEXP, SEXP segment_lengthSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type phylogeny(phylogenySEXP);
    Rcpp::traits::input_parameter< List >::type mode_phylogeny(mode_phylogenySEXP);
    Rcpp::traits::input_parameter< List >::type models(modelsSEXP);
    Rcpp::traits::input_parameter< List >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type start_mode(start_modeSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type segment_length(segment_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_over_interval_phylogeny(phylogeny, mode_phylogeny, models, sequence, start_mode, rate, segment_length, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// simulate_with_shared_substitution_heterogeneity
List simulate_with_shared_substitution_heterogeneity(List phylogeny, List switching_model, List substitution_models, List sequence, unsigned long long start_mode, double rate, double switching_rate, double segment_length, double tolerance);
RcppExport SEXP _PalantiR_simulate_with_shared_substitution_heterogeneity(SEXP phylogenySEXP, SEXP switching_modelSEXP, SEXP substitution_modelsSEXP, SEXP sequenceSEXP, SEXP start_modeSEXP, SEXP rateSEXP, SEXP switching_rateSEXP, SEXP segment_lengthSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type phylogeny(phylogenySEXP);
    Rcpp::traits::input_parameter< List >::type switching_model(switching_modelSEXP);
    Rcpp::traits::input_parameter< List >::type substitution_models(substitution_modelsSEXP);
    Rcpp::traits::input_parameter< List >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type start_mode(start_modeSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type switching_rate(switching_rateSEXP);
    Rcpp::traits::input_parameter< double >::type segment_length(segment_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_with_shared_substitution_heterogeneity(phylogeny, switching_model, substitution_models, sequence, start_mode, rate, switching_rate, segment_length, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// simulate_with_shared_time_heterogeneity
List simulate_with_shared_time_heterogeneity(List phylogeny, List switching_model, List substitution_models, List sequence, unsigned long long start_mode, double rate, double switching_rate, double segment_length, double tolerance);
RcppExport SEXP _PalantiR_simulate_with_shared_time_heterogeneity(SEXP phylogenySEXP, SEXP switching_modelSEXP, SEXP substitution_modelsSEXP, SEXP sequenceSEXP, SEXP start_modeSEXP, SEXP rateSEXP, SEXP switching_rateSEXP, SEXP segment_lengthSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type phylogeny(phylogenySEXP);
    Rcpp::traits::input_parameter< List >::type switching_model(switching_modelSEXP);
    Rcpp::traits::input_parameter< List >::type substitution_models(substitution_modelsSEXP);
    Rcpp::traits::input_parameter< List >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type start_mode(start_modeSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type switching_rate(switching_rateSEXP);
    Rcpp::traits::input_parameter< double >::type segment_length(segment_lengthSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_with_shared_time_heterogeneity(phylogeny, switching_model, substitution_models, sequence, start_mode, rate, switching_rate, segment_length, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// compare_modes
bool compare_modes(const List& substitution_models, const List& intervals);
RcppExport SEXP _PalantiR_compare_modes(SEXP substitution_modelsSEXP, SEXP intervalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type substitution_models(substitution_modelsSEXP);
    Rcpp::traits::input_parameter< const List& >::type intervals(intervalsSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_modes(substitution_models, intervals));
    return rcpp_result_gen;
END_RCPP
}
// Sequence
List Sequence(std::string sequence, std::string type, unsigned long long mode);
RcppExport SEXP _PalantiR_Sequence(SEXP sequenceSEXP, SEXP typeSEXP, SEXP modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type mode(modeSEXP);
    rcpp_result_gen = Rcpp::wrap(Sequence(sequence, type, mode));
    return rcpp_result_gen;
END_RCPP
}
// as_compound
List as_compound(List sequence, unsigned long long mode);
RcppExport SEXP _PalantiR_as_compound(SEXP sequenceSEXP, SEXP modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type mode(modeSEXP);
    rcpp_result_gen = Rcpp::wrap(as_compound(sequence, mode));
    return rcpp_result_gen;
END_RCPP
}
// sample_sequence
List sample_sequence(List model, unsigned long long length);
RcppExport SEXP _PalantiR_sample_sequence(SEXP modelSEXP, SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type model(modelSEXP);
    Rcpp::traits::input_parameter< unsigned long long >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_sequence(model, length));
    return rcpp_result_gen;
END_RCPP
}
// as_amino_acid
std::vector<std::string> as_amino_acid(std::vector<std::string>& codons);
RcppExport SEXP _PalantiR_as_amino_acid(SEXP codonsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type codons(codonsSEXP);
    rcpp_result_gen = Rcpp::wrap(as_amino_acid(codons));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PalantiR_GeneticCode", (DL_FUNC) &_PalantiR_GeneticCode, 0},
    {"_PalantiR_HasegawaKishinoYano", (DL_FUNC) &_PalantiR_HasegawaKishinoYano, 3},
    {"_PalantiR_GeneralTimeReversible", (DL_FUNC) &_PalantiR_GeneralTimeReversible, 2},
    {"_PalantiR_MutationSelection", (DL_FUNC) &_PalantiR_MutationSelection, 5},
    {"_PalantiR_CoEvolution", (DL_FUNC) &_PalantiR_CoEvolution, 7},
    {"_PalantiR_MarkovModulatedMutationSelection", (DL_FUNC) &_PalantiR_MarkovModulatedMutationSelection, 2},
    {"_PalantiR_Phylogeny", (DL_FUNC) &_PalantiR_Phylogeny, 2},
    {"_PalantiR_equilibrium_to_fitness", (DL_FUNC) &_PalantiR_equilibrium_to_fitness, 2},
    {"_PalantiR_simulate_over_phylogeny", (DL_FUNC) &_PalantiR_simulate_over_phylogeny, 4},
    {"_PalantiR_phylogeny_to_intervals", (DL_FUNC) &_PalantiR_phylogeny_to_intervals, 2},
    {"_PalantiR_simulate_over_interval_phylogeny", (DL_FUNC) &_PalantiR_simulate_over_interval_phylogeny, 8},
    {"_PalantiR_simulate_with_shared_substitution_heterogeneity", (DL_FUNC) &_PalantiR_simulate_with_shared_substitution_heterogeneity, 9},
    {"_PalantiR_simulate_with_shared_time_heterogeneity", (DL_FUNC) &_PalantiR_simulate_with_shared_time_heterogeneity, 9},
    {"_PalantiR_compare_modes", (DL_FUNC) &_PalantiR_compare_modes, 2},
    {"_PalantiR_Sequence", (DL_FUNC) &_PalantiR_Sequence, 3},
    {"_PalantiR_as_compound", (DL_FUNC) &_PalantiR_as_compound, 2},
    {"_PalantiR_sample_sequence", (DL_FUNC) &_PalantiR_sample_sequence, 2},
    {"_PalantiR_as_amino_acid", (DL_FUNC) &_PalantiR_as_amino_acid, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_PalantiR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
