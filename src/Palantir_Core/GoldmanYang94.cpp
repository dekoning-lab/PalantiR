#include "GoldmanYang94.hpp"

mat Palantir::GoldmanYang94::transition(
        const vec& equilibrium,
        double omega,
        double kappa,
        const GeneticCode& g)
{
    if(equilibrium.n_elem != g.size) {
        throw invalid_argument("GoldmanYang94 requires one equilibrium "
                               "frequency per sense codon");
    }

    mat codon_transition(g.size, g.size, fill::zeros);
    for(const Codon& i : g) {
        for(const Codon& j : g) {
            if(i == j || Codon::_distance(i, j) != 1) {
                continue;
            }

            pair<ullong, ullong> change = Codon::_substitution(i, j);
            double relative_rate = equilibrium[j.index];
            if(Nucleotide::transition(change.first, change.second)) {
                relative_rate *= kappa;
            }
            if(!Codon::_synonymous(i, j)) {
                relative_rate *= omega;
            }
            codon_transition.at(i.index, j.index) = relative_rate;
        }
    }

    codon_transition.diag() = -sum(codon_transition, 1);
    return codon_transition;
}
