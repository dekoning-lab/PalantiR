#include <RcppArmadillo.h>
using namespace Rcpp;

#include "RcppUtil.hpp"
#include "GeneticCode.hpp"

void decorate_substitutions(List& substitutions, std::string type);
