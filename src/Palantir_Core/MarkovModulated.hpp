#ifndef MarkovModulated_hpp
#define MarkovModulated_hpp

#include "Util.hpp"

#include "GeneticCode.hpp"

namespace Palantir
{
    namespace MarkovModulated
    {
        // FIX (2026-08-20): the equilibrium() declaration was removed along
        // with its (dead, incorrect) definition; use
        // MarkovModel::solve_equilibrium on the transition matrix instead.

        mat transition(
                const vector<mat>& substitution_transition,
                const mat& exchangeability,
                const vec& switching_equilibrium,
                const GeneticCode& g);

    }
}
#endif /* MarkovModulated_hpp */
