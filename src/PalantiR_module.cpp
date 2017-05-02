#include <RcppArmadillo.h>
#include "Phylogeny.hpp"
using namespace Rcpp;


RCPP_EXPOSED_CLASS(Phylogeny)

RCPP_MODULE(PalantiR_module)
{
    class_<Palantir::Phylogeny>("PhylogenyClass")
    .constructor<std::string>()
    .field_readonly("n_nodes", &Palantir::Phylogeny::n_nodes)
    .method("to_JSON", &Palantir::Phylogeny::to_JSON);
}
