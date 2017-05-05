#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

std::string get_genetic_code_name();

bool has_class(List a, std::string cl);

CharacterVector predicate(arma::uvec predicate, string yes, string no);
