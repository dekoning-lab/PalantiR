#ifndef CodonModel_hpp
#define CodonModel_hpp

#include "Util.hpp"
#include "GeneticCode.hpp"

namespace Palantir
{
    namespace CodonModel
    {
        // PalantiR historically calls the conventional all-substitution scale
        // "substitution".  "standard" is accepted as a user-facing alias.
        string canonical_scaling_type(const string& scaling_type);

        // Per-state outflux in the substitution class used to denominate branch
        // lengths.  Q is expected to be a single-codon generator in the active
        // genetic code's sense-codon order.
        vec class_outflux(
                const mat& transition,
                const string& scaling_type,
                const GeneticCode& g);

        // Stationary expected event rate for the set of transitions that
        // determines branch-length scaling.
        double scaling(
                const vec& equilibrium,
                const mat& transition,
                const string& scaling_type,
                const GeneticCode& g);
    }
}

#endif /* CodonModel_hpp */
