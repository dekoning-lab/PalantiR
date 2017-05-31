#ifndef Simulate_hpp
#define Simulate_hpp

#include "Palantir.hpp"
#include "Phylogeny.hpp"
#include "GeneticCode.hpp"
#include "SubstitutionHistory.hpp"
#include "IntervalHistory.hpp"
#include "SiteSimulation.hpp"

// Only scaling:
#include "MutationSelection.hpp"

namespace Palantir
{
    namespace Simulate
    {

        vector<double> poisson(
                double time = 1,
                double rate = 1);

        vector<ullong> steps(
                const mat& sampling,
                ullong start = 0,
                ullong steps = 1);

        SubstitutionHistory over_time(
                const mat& transition,
                const mat& sampling,
                ullong start = 0,
                double time = 1,
                double rate = 1);

        SiteSimulation over_phylogeny(
                const Phylogeny& tree,
                const mat& transition,
                const mat& sampling,
                ullong start = 0,
                double rate = 1);

        vector<SiteSimulation> sequence_over_phylogeny(
                const Phylogeny& tree,
                const mat& transition,
                const mat& sampling,
                const uvec& sequence,
                double rate = 1);

        vector<SiteSimulation> sequence_over_intervals(
                const Phylogeny& tree,
                const vector<IntervalHistory>& tree_intervals,
                const vector<vec>& equilibrium,
                const vector<mat>& transition,
                const vector<mat>& sampling,
                const uvec& sequence,
                ullong start_mode,
                const GeneticCode& g,
                double rate = 1,
                double segment_length = 0.001,
                double tolerance = 0.001,
                string scaling_type = "synonymous");

        vector<IntervalHistory> switching_intervals(
                const Phylogeny& tree,
                const mat& switching_transition,
                const mat& switching_sampling,
                ullong start_mode,
                double rate = 1);

        vector<vector<IntervalHistory>> switching_poisson(
                const Phylogeny& tree,
                const mat& switching_sampling,
                ullong n_sites,
                ullong start_mode,
                double rate = 1);
    }
}

#endif /* Simulate_hpp */
