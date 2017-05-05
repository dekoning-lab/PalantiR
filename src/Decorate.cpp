#include "Decorate.hpp"

void decorate_substitutions(List& substitutions, std::string type)
{
    Palantir::GeneticCode g(get_genetic_code_name());
    arma::uvec state_from = substitutions["from"];
    arma::uvec state_to = substitutions["to"];

    if (type == "codon") {
        substitutions["codon_from"] = Palantir::Codon::to_string(state_from, g);
        substitutions["codon_to"] = Palantir::Codon::to_string(state_to, g);
        substitutions["amino_acid_from"] = Palantir::Codon::to_amino_acid(state_from, g);
        substitutions["amino_acid_to"] = Palantir::Codon::to_amino_acid(state_to, g);

        arma::uvec synonymous = Palantir::Codon::synonymous(state_from, state_to, g);
        substitutions["synonymous"] = LogicalVector(synonymous.begin(), synonymous.end());
        substitutions["color"] = predicate(synonymous, "#16A35E", "#E84B2A");
    } else {
        stop("Unknown substitution history type");
    }
}
