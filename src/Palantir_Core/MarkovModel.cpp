#include "MarkovModel.hpp"

mat Palantir::MarkovModel::random_transition(const ullong size)
{
    mat transition_rates(size, size, fill::randu);
    transition_rates.diag() = zeros(size);
    
    vec row_sums = sum(transition_rates, 1);
    transition_rates.diag() = -row_sums;
    
    return transition_rates;
}

vec Palantir::MarkovModel::solve_equilibrium(const mat& transition)
{
    ullong n = transition.n_cols;
    mat U;
    vec s;
    mat V;
    svd(U, s, V, transition);
    
    vec pi = U.col(n-1);
    
    pi /= sum(pi);
    
    return pi;
}

