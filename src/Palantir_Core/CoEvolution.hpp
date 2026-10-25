#ifndef CoEvolution_hpp
#define CoEvolution_hpp

#include "Util.hpp"

#include "GeneticCode.hpp"
#include "MutationSelection.hpp"

#include <armadillo>


namespace Palantir{
    namespace CoEvolution
    {
        vec equilibrium(
                ullong population_size,
                double mutation_rate,
                const vec& nucleotide_equilibrium,
                const vec& fitness_1,
                const vec& fitness_2,
                const mat& delta,
                const GeneticCode& g,
                FixationModel model = FixationModel::Diffusion);

        mat transition(
                ullong population_size,
                double mutation_rate,
                const mat& nucleotide_transition,
                const vec& fitness_1,
                const vec& fitness_2,
                const mat& delta,
                const GeneticCode& g,
                FixationModel model = FixationModel::Diffusion);

        double scaling(
                const vec& equilibrium,
                const mat& transition,
                const string scaling_type,
                const GeneticCode& g);
    }
}
#endif /* CoEvolution_hpp */
