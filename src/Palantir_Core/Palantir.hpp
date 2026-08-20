#ifndef Palantir_hpp
#define Palantir_hpp

#include "Util.hpp"

    bool approximately_equal(const vec& a, const vec& b, double tolerance);
    double rmsd(const vec& a, const vec& b);
    double rnd_exp(double rate);

namespace Palantir {
    // Declared in GeneticCode.hpp; only referenced here.
    class GeneticCode;

    // FIX (2026-08-20): the engine used to be declared here as
    //     static random_device rd;
    //     static mt19937 rng(rd());
    // i.e. at namespace scope in a header, so every translation unit got its
    // OWN generator, each seeded independently from random_device. Waiting
    // times (Simulate.cpp) and jump/state draws (Palantir.cpp) came from
    // different streams, nothing could be seeded, and no simulation was
    // reproducible. There is now exactly one engine, defined once in
    // Palantir_Core/Palantir.cpp and reseedable via seed_rng() (exposed to R as
    // set_palantir_seed()). Call sites are unchanged: `Palantir::rng` still
    // names the engine and rnd_exp(rate) still works.
    extern mt19937 rng;

    void seed_rng(unsigned int seed);

// Fallback for approx_equal(a, b, "absdiff")

    mat sampling(const mat& transition_rates);

    uvec sample_sequence(const vec& equilibrium, ullong length);

    vec equilibrium_to_fitness(const vec& equilibrium, double population_size, double epsilon = 1e-316);

    // FIX (2026-08-20, M2): degeneracy-aware inversion. The two-argument form
    // above inverts pi -> psi ignoring the mutational weight of each amino
    // acid's codons, so the amino-acid marginal of the resulting
    // mutation-selection equilibrium is the codon-degeneracy-tilted transform
    // of the target rather than the target itself. Given the nucleotide
    // equilibrium of the mutation model, this overload divides out that weight
    // so that the realised marginal equals `equilibrium` exactly. The
    // two-argument form is left untouched for backward compatibility.
    vec amino_acid_mutational_weight(const vec& nucleotide_equilibrium, const GeneticCode& g);

    vec equilibrium_to_fitness(
            const vec& equilibrium,
            double population_size,
            const vec& nucleotide_equilibrium,
            const GeneticCode& g,
            double epsilon = 1e-316);
}


#endif /* Palantir_hpp */
