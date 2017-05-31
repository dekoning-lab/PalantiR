#ifndef MarkovModel_hpp
#define MarkovModel_hpp

#include "Util.hpp"

namespace Palantir
{
    namespace MarkovModel 
    {
        vec solve_equilibrium(const mat& transition);
        mat random_transition(const ullong size);
    }
}

#endif /* MarkovModel_hpp */
