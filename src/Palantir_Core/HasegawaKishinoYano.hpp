#ifndef HasegawaKishinoYano_hpp
#define HasegawaKishinoYano_hpp

#include "Util.hpp"
#include "Palantir.hpp"
#include "GeneticCode.hpp"

namespace Palantir 
{
    namespace HasegawaKishinoYano
    {
        mat transition(
                const vec& equilibrium,
                double transition_rate = 1,
                double transversion_rate = 1);
    }
}

#endif /* HasegawaKishinoYano_hpp */
