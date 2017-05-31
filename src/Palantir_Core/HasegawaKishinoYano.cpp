#include "HasegawaKishinoYano.hpp"

mat Palantir::HasegawaKishinoYano::transition(
        const vec& equilibrium,
        double transversion_rate,
        double transition_rate)
{
    ullong size = equilibrium.n_elem;
    if(!(size == 4)) {
        throw logic_error("Equilibrium distribution of a nucleotide model should have 4 states");
    }
    
    mat transition(size, size, fill::zeros);

    for(ullong i = 0; i < Nucleotide::size; i++) {
        for(ullong j = 0; j < Nucleotide::size; j++) {
            if (i != j) {
                double r = transversion_rate;
                if(Nucleotide::transition(i, j)) {
                    r = transition_rate;
                }
                transition.at(i, j) = r * equilibrium[j];
            }
        }
    }
    
    vec row_sums = sum(transition, 1);
    transition.diag() = -row_sums;
    
    return transition;
}
