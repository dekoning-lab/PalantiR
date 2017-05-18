#include "RcppPalantir.hpp"

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
    } else if(type == "compound_codon") {
        uvec codon_from  = Palantir::CompoundCodon::to_codon_index(state_from, g);
        uvec codon_to = Palantir::CompoundCodon::to_codon_index(state_to, g);
        uvec model_from = Palantir::CompoundCodon::to_mode_index(state_from, g);
        uvec model_to = Palantir::CompoundCodon::to_mode_index(state_to, g);

        substitutions["model_from"] = model_from;
        substitutions["model_to"] = model_to;
        substitutions["codon_from"] = Palantir::Codon::to_string(codon_from, g);
        substitutions["codon_to"] = Palantir::Codon::to_string(codon_to, g);

        arma::uvec synonymous = Palantir::Codon::synonymous(codon_from, codon_to, g);
        substitutions["synonymous"] = LogicalVector(synonymous.begin(), synonymous.end());

        arma::uvec model_switch(model_from.n_elem);
        for(unsigned long long i = 0; i < model_from.n_elem; i++) {
            model_switch[i] = (model_to[i] != model_from[i]);
        }
        substitutions["model_switch"] = LogicalVector(model_switch.begin(), model_switch.end());

        CharacterVector color = predicate(synonymous, "#16A35E", "#E84B2A");
        for(unsigned long long i = 0; i < model_switch.n_elem; i++) {
            if(model_switch[i]) {
                color[i] = "#002FCF";
            }
        }
        substitutions["color"] = color;
    } else if(type == "codon_pair") {
        substitutions["pair_from"] = Palantir::CodonPair::to_string(state_from, g);
        substitutions["pair_to"] = Palantir::CodonPair::to_string(state_to, g);

        uvec first_from = Palantir::CodonPair::first(state_from, g);
        uvec second_from = Palantir::CodonPair::second(state_from, g);

        uvec first_to = Palantir::CodonPair::first(state_to, g);
        uvec second_to = Palantir::CodonPair::second(state_to, g);

        vector<string> first_aa_from = Palantir::Codon::to_amino_acid(first_from, g);
        vector<string> second_aa_from = Palantir::Codon::to_amino_acid(second_from, g);

        vector<string> first_aa_to = Palantir::Codon::to_amino_acid(first_to, g);
        vector<string> second_aa_to = Palantir::Codon::to_amino_acid(second_to, g);

        vector<string> pair_aa_from(first_aa_from.size());
        vector<string> pair_aa_to(first_aa_to.size());
        for(ullong i = 0; i < first_aa_from.size(); i++) {
            pair_aa_from[i] = first_aa_from[i] + "," + second_aa_from[i];
            pair_aa_to[i] = first_aa_to[i] + "," + second_aa_to[i];
        }
        substitutions["pair_amino_acid_from"] = pair_aa_from;
        substitutions["pair_amino_acid_to"] = pair_aa_to;

        uvec first_synonymous = Palantir::Codon::synonymous(first_from, first_to, g);
        substitutions["first_synonymous"] = LogicalVector(first_synonymous.begin(), first_synonymous.end());

        uvec second_synonymous = Palantir::Codon::synonymous(second_from, second_to, g);
        substitutions["second_synonymous"] = LogicalVector(second_synonymous.begin(), second_synonymous.end());

        substitutions["first_color"] = predicate(first_synonymous, "#16A35E", "#E84B2A");
        substitutions["second_color"] = predicate(second_synonymous, "#16A35E", "#E84B2A");
    } else {
        stop("Unknown substitution history type");
    }
}
