#ifndef CodonModel_hpp
#define CodonModel_hpp

#include "Util.hpp"
#include "GeneticCode.hpp"

namespace Palantir
{
    namespace CodonModel
    {
        // "synonymous" is the user-facing dS gauge; "dS" and "ds" are
        // explicit aliases. The historical synonymous-events-per-codon
        // currency remains available under the unambiguous name
        // "synonymous-per-codon".
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

        // Denominator that converts a neutral single-codon generator to the
        // dS gauge: Q_dS = Q_neutral / neutral_dS_scaling(...). The resulting
        // neutral total rate is three substitutions per codon.
        double neutral_dS_scaling(
                const vec& neutral_equilibrium,
                const mat& neutral_transition,
                const GeneticCode& g);

        // Per-state rate at which the neutral reference consumes dS time.
        // `neutral_transition` must already be expressed in the model's dS
        // gauge, so its stationary mean total rate is three per codon. Dividing
        // the total outflux by three converts it to a per-nucleotide clock.
        vec neutral_dS_outflux(
                const mat& neutral_transition,
                const GeneticCode& g);
    }
}

#endif /* CodonModel_hpp */
