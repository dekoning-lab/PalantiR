#ifndef MutationSelection_hpp
#define MutationSelection_hpp

#include "Util.hpp"
#include "Palantir.hpp"

#include "GeneticCode.hpp"
#include "HasegawaKishinoYano.hpp"


namespace Palantir 
{
    namespace MutationSelection
    {
        vec equilibrium(
                ullong population_size,
                double mutation_rate,
                const vec& nucleotide_equilibrium,
                const vec& fitness,
                const GeneticCode& g,
                FixationModel model = FixationModel::Diffusion);

        mat transition(
                ullong population_size,
                double mutation_rate,
                const mat& nucleotide_transition,
                const vec& fitness,
                const GeneticCode& g,
                FixationModel model = FixationModel::Diffusion);

        double scaling(
                const vec& equilibrium,
                const mat& transition,
                const string scaling_type,
                const GeneticCode& g);

        double fixation_probability(
                ullong population_size,
                double selection,
                FixationModel model = FixationModel::Diffusion);
    }
}

#endif /* MutationSelection_hpp */
