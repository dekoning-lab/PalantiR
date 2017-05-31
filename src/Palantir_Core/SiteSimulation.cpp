#include "SiteSimulation.hpp"

Palantir::SiteSimulation::SiteSimulation(const vector<SubstitutionHistory>& substitutions, const uvec& states):
    substitutions(substitutions),
    states(states)
{
}

uvec Palantir::SiteSimulation::leaf_states(const Phylogeny& p) const
{
    vector<reference_wrapper<const Phylogeny::Node> > leaves = p.leaf_traversal();
    uvec leaf_states(leaves.size());

    for(ullong i = 0; i < leaves.size(); i++) {
        const Phylogeny::Node& leaf = leaves[i].get();
        leaf_states[i] = states[leaf.index];
    }
    return leaf_states;
}
