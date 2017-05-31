#include "GeneralTimeReversible.hpp"

mat Palantir::GeneralTimeReversible::transition(
        const vec& equilibrium,
        const mat& exchangeability,
        double rate)
{
    ullong size = equilibrium.n_elem;

    mat diff = abs(trans(exchangeability) - exchangeability);
    bool symmetric = all(vectorise(diff) <= 0);
    if(!symmetric) {
        throw logic_error("Exchangeability matrix must be symmetric");
    }
    
    mat transition(size, size, fill::zeros);

    for(ullong i = 0; i < size; i++) {
        for(ullong j = 0; j < size; j++) {
            if (i != j) {
                transition(i,j) = equilibrium[j] * exchangeability.at(i, j) * rate;
            }
        }
    }
    
    vec row_sums = sum(transition, 1);
    transition.diag() = -row_sums;
    
    return transition;
}

