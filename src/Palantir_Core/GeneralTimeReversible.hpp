#ifndef GeneralTimeReversible_hpp
#define GeneralTimeReversible_hpp

#include "Util.hpp"
#include "Palantir.hpp"

using namespace std;

namespace Palantir
{
    namespace GeneralTimeReversible 
    {
        mat transition(
                const vec& equilibrium,
                const mat& exchangeability,
                double rate = 1.0);
    }
}

#endif /* GeneralTimeReversible_hpp */
