#include "MarkovModulated.hpp"

// FIX (2026-08-20): MarkovModulated::equilibrium was deleted. It was dead code
// (Models.cpp uses MarkovModel::solve_equilibrium) and mathematically wrong:
// w_i * pi_i(c) is not the stationary distribution of the modulated chain when
// the component equilibria differ, which is the whole point of the model.

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

    // FIX (2026-08-20): the requested mode occupancy was ignored -- the
    // off-diagonal (switching) blocks used the bare exchangeability, whose
    // symmetry makes the mode marginal uniform whatever occupancy was asked
    // for -- while switching_equilibrium instead DIVIDED the diagonal
    // substitution blocks, inflating every mode's substitution rates by
    // exactly 1/w_i (factor 1.25 at w = 0.8), so branch lengths meant nothing.
    // Build the mode process as a proper GTR instead: the switching rate from
    // mode i to mode j is S_ij * w_j (S the symmetric exchangeability, w the
    // requested occupancy), which is reversible with stationary distribution w,
    // and leave each mode's substitution block as the component model's matrix
    // unchanged, preserving the components' rate normalisation.
    for (ullong i = 0; i < n_modes; i++) {
        for (ullong j = 0; j < n_modes; j++) {
            if (i != j) {
                transition.submat(
                        i * size,
                        j * size,
                        ((i + 1) * size) - 1,
                        ((j + 1) * size) - 1) = eye(size, size) *
                            (exchangeability.at(i, j) * switching_equilibrium[j]);
            } else {
                transition.submat(
                        i * size,
                        i * size,
                        ((i + 1) * size) - 1,
                        ((i + 1) * size) - 1) = substitution_transition_rates[i];
            }

        }
    }

    transition.diag() = zeros(state_space);
    vec row_sums = sum(transition, 1);
    transition.diag() = -row_sums;
    
    return transition;
}
