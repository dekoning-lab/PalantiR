// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Sequence
List Sequence(std::string sequence, std::string type, unsigned long long mode);
RcppExport SEXP PalantiR_Sequence(SEXP sequenceSEXP, SEXP typeSEXP, SEXP modeSEXP) {
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
// HasegawaKishinoYano
List HasegawaKishinoYano(arma::vec equilibrium, double transition_rate, double transversion_rate);
RcppExport SEXP PalantiR_HasegawaKishinoYano(SEXP equilibriumSEXP, SEXP transition_rateSEXP, SEXP transversion_rateSEXP) {
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
// MutationSelection
List MutationSelection(unsigned long long population_size, double mutation_rate, arma::vec nucleotide_equilibrium, arma::mat nucleotide_transition, arma::vec fitness, std::string scaling_type);
RcppExport SEXP PalantiR_MutationSelection(SEXP population_sizeSEXP, SEXP mutation_rateSEXP, SEXP nucleotide_equilibriumSEXP, SEXP nucleotide_transitionSEXP, SEXP fitnessSEXP, SEXP scaling_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned long long >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nucleotide_equilibrium(nucleotide_equilibriumSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type nucleotide_transition(nucleotide_transitionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fitness(fitnessSEXP);
    Rcpp::traits::input_parameter< std::string >::type scaling_type(scaling_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(MutationSelection(population_size, mutation_rate, nucleotide_equilibrium, nucleotide_transition, fitness, scaling_type));
    return rcpp_result_gen;
END_RCPP
}
// decorate_codon_substitutions
DataFrame decorate_codon_substitutions(DataFrame substitutions);
RcppExport SEXP PalantiR_decorate_codon_substitutions(SEXP substitutionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type substitutions(substitutionsSEXP);
    rcpp_result_gen = Rcpp::wrap(decorate_codon_substitutions(substitutions));
    return rcpp_result_gen;
END_RCPP
}
// Phylogeny
List Phylogeny(std::string newick);
RcppExport SEXP PalantiR_Phylogeny(SEXP newickSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type newick(newickSEXP);
    rcpp_result_gen = Rcpp::wrap(Phylogeny(newick));
    return rcpp_result_gen;
END_RCPP
}
// equilibrium_to_fitness
arma::vec equilibrium_to_fitness(arma::vec equilibrium, unsigned long long population_size);
RcppExport SEXP PalantiR_equilibrium_to_fitness(SEXP equilibriumSEXP, SEXP population_sizeSEXP) {
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
List simulate_over_phylogeny(List tree, List substitution_model, arma::uvec sequence, double rate);
RcppExport SEXP PalantiR_simulate_over_phylogeny(SEXP treeSEXP, SEXP substitution_modelSEXP, SEXP sequenceSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< List >::type substitution_model(substitution_modelSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_over_phylogeny(tree, substitution_model, sequence, rate));
    return rcpp_result_gen;
END_RCPP
}
