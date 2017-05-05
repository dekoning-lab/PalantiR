#include <RcppArmadillo.h>
using namespace Rcpp;

#include "Palantir_Core/GeneticCode.hpp"
#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/Phylogeny.hpp"
#include "Palantir_Core/SiteSimulation.hpp"

std::string get_genetic_code_name();

bool has_class(List a, std::string cl);

CharacterVector predicate(arma::uvec predicate, string yes, string no);

void decorate_substitutions(List& substitutions, std::string type);

CharacterMatrix get_alignment(const vector<Palantir::SiteSimulation>& sims, const Palantir::Phylogeny& p, std::string type);
