#ifndef MarkovModulated_hpp
#define MarkovModulated_hpp

#include "Util.hpp"

#include "GeneticCode.hpp"

namespace Palantir
{
    namespace MarkovModulated
    {
        vec equilibrium(
                const vec& switching_equilibrium,
                const vector<vec>& substitution_equilibrium,
                const GeneticCode& g);

        mat transition(
                const vector<mat>& substitution_transition,
                const mat& exchangeability,
                const vec& switching_equilibrium,
                const GeneticCode& g);

    }
}
#endif /* MarkovModulated_hpp */
