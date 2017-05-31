#ifndef Palantir_hpp
#define Palantir_hpp

#include "Util.hpp"

    bool approximately_equal(const vec& a, const vec& b, double tolerance);
    double rmsd(const vec& a, const vec& b);
    double rnd_exp(double rate);

namespace Palantir {
    static random_device rd;
    static mt19937 rng(rd());

// Fallback for approx_equal(a, b, "absdiff")

    mat sampling(const mat& transition_rates);

    uvec sample_sequence(const vec& equilibrium, ullong length);

    vec equilibrium_to_fitness(const vec& equilibrium, ullong population_size, double epsilon = 1e-316);
}


#endif /* Palantir_hpp */
