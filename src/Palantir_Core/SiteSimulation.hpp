#ifndef SiteSimulation_hpp
#define SiteSimulation_hpp

#include "Util.hpp"
#include "SubstitutionHistory.hpp"
#include "Phylogeny.hpp"

namespace Palantir
{
    class SiteSimulation 
    {
        public:

            const vector<SubstitutionHistory> substitutions;
            const uvec states;

            SiteSimulation(const vector<SubstitutionHistory>& substitutions, const uvec& states);

            uvec leaf_states(const Phylogeny& p) const;
    };
}

#endif //SiteSimulation_hpp
