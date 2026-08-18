#include "HasegawaKishinoYano.hpp"

// FIX (2026-08-17): the parameter order here disagreed with the declaration in
// HasegawaKishinoYano.hpp, which is (equilibrium, transition_rate,
// transversion_rate). C++ binds positionally, so every caller asking for
// transition_rate = k silently got transversion_rate = k. Reordered to match the
// header. All other call sites use defaults only, so nothing else changes.
mat Palantir::HasegawaKishinoYano::transition(
        const vec& equilibrium,
        double transition_rate,
        double transversion_rate)
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
