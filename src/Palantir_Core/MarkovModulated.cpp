#include "MarkovModulated.hpp"

vec Palantir::MarkovModulated::equilibrium(
        const vec& switching_equilibrium,
        const vector<vec>& substitution_equilibrium,
        const GeneticCode& g)
{
    ullong n_modes = switching_equilibrium.n_elem;
    
    vec equilibrium_distribution(n_modes * g.size, fill::zeros);
    
    for(ullong i = 0; i < n_modes; i++) {
        const vec mode_equilibrium(substitution_equilibrium[i]);
        for(const Codon c : g) {
            ullong j = c.index;
            equilibrium_distribution((i * g.size) + j) = mode_equilibrium[j] * switching_equilibrium[i];
        }
    }
    return equilibrium_distribution;
}

mat Palantir::MarkovModulated::transition(
        const vector<mat>& substitution_transition_rates,
        const mat& exchangeability,
        const vec& switching_equilibrium,
        const GeneticCode& g)
{
    ullong n_modes = switching_equilibrium.n_elem;
    ullong size = g.size;
    ullong state_space = n_modes * size;
    
    mat transition(state_space, state_space, fill::zeros);

    for (ullong i = 0; i < n_modes; i++) {
        for (ullong j = 0; j < n_modes; j++) {
            if (i != j) {
                transition.submat(
                        i * size,
                        j * size,
                        ((i + 1) * size) - 1,
                        ((j + 1) * size) - 1) = eye(size, size) * exchangeability.at(i, j);
            } else {
                transition.submat(
                        i * size,
                        i * size,
                        ((i + 1) * size) - 1,
                        ((i + 1) * size) - 1) = (substitution_transition_rates[i] / switching_equilibrium[i]);
            }
            
        }
    }

    transition.diag() = zeros(state_space);
    vec row_sums = sum(transition, 1);
    transition.diag() = -row_sums;
    
    return transition;
}
