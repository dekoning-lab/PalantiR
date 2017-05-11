#include <RcppArmadillo.h>
#include <set>
using namespace Rcpp;

#include "Palantir_Core/GeneticCode.hpp"
#include "Palantir_Core/Palantir.hpp"
#include "Palantir_Core/Phylogeny.hpp"
#include "Palantir_Core/SiteSimulation.hpp"
#include "Palantir_Core/Util.hpp"

typedef unsigned long long ullong;

std::string get_genetic_code_name();

std::string get_attr(List l, std::string a);

bool has_class(List a, std::string cl);

CharacterVector predicate(arma::uvec predicate, string yes, string no);

void decorate_substitutions(List& substitutions, std::string type);

CharacterMatrix get_alignment(const vector<Palantir::SiteSimulation>& sims, const Palantir::Phylogeny& p, std::string type);

List interval_histories_to_list(const std::vector<Palantir::IntervalHistory>& tree_intervals, const Palantir::Phylogeny& p);
List site_simulations_to_list(const std::vector<Palantir::SiteSimulation>& sims, const Palantir::Phylogeny& p);

bool compare_modes(const List& substitution_models, const List& intervals);
