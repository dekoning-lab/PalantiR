#ifndef GoldmanYang94_hpp
#define GoldmanYang94_hpp

#include "Util.hpp"
#include "GeneticCode.hpp"

namespace Palantir
{
    namespace GoldmanYang94
    {
        // Goldman-Yang (1994) codon generator.  `equilibrium` is the stationary
        // distribution over the active code's sense codons.  Only single-
        // nucleotide changes are permitted; kappa multiplies transitions and
        // omega multiplies non-synonymous changes.
        mat transition(
                const vec& equilibrium,
                double omega,
                double kappa,
                const GeneticCode& g);
    }
}

#endif /* GoldmanYang94_hpp */
