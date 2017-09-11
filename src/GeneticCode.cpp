#include <RcppArmadillo.h>

#include "Palantir_Core/GeneticCode.hpp"

#include "RcppPalantir.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
StringVector GeneticCode() {
    Palantir::GeneticCode g(get_genetic_code_name());

    StringVector gc(g.size);
    std::vector<std::string> codons(g.size);
    for (ullong i = 0; i < g.size; i++) {
        ullong aa = g.codons[i].amino_acid;
        gc[i] = Palantir::AminoAcid::AminoAcids[aa];

        codons[i] = g.sequence[i];
    }

    gc.names() = codons;
    return gc;
}
