#include "test.hpp"

TEST_CASE("PopulationSwitching")
{
    vec nucleotide_equilibrium("0.25 0.25 0.25 0.25");
    vec switching_equilibrium("0.25 0.25 0.25 0.25");
    vec fitness(test_fitness_1);
    GeneticCode g("Standard nuclear");

    Phylogeny tree(test_tree);
    mat exchangeability("0 0.2 0.1 0.3; 0.2 0 0.1 0.2; 0.1 0.1 0 0.1; 0.3 0.2 0.1 0");
    vec population_sizes("10 100 1000 10000");

    mat nucleotide_transition = HasegawaKishinoYano::transition(nucleotide_equilibrium);
    mat switching_transition = GeneralTimeReversible::transition(switching_equilibrium, exchangeability);

    ullong n_modes = population_sizes.n_elem;

    ullong state_space = g.size * n_modes;
    vector<vec> substitution_equilibrium;
    vector<mat> substitution_transition;
    vector<double> substitution_rho;

    for(ullong i = 0; i < n_modes; i++) {
        vec codon_equilibrium = MutationSelection::equilibrium(
                population_sizes.at(i), 1e-8, nucleotide_equilibrium, fitness, g);
        mat codon_transition = MutationSelection::transition(
                population_sizes.at(i), 1e-8, nucleotide_transition, fitness, g);
        double codon_rho = MutationSelection::scaling(
                codon_equilibrium, codon_transition, "synonymous", g);
        codon_transition /= codon_rho;

        substitution_equilibrium.push_back(codon_equilibrium);
        substitution_transition.push_back(codon_transition);
        substitution_rho.push_back(codon_rho);
    }

    mat transition = MarkovModulated::transition(
            substitution_transition, exchangeability, switching_equilibrium, g);

    // This works rather well
    vec equilibrium = MarkovModel::solve_equilibrium(transition);

    mat sampling = Palantir::sampling(transition);

    SECTION("Switching Model check") {
        REQUIRE(sum(equilibrium) - 1.0 < SMALL);


        REQUIRE(all(vectorise(abs(sum(transition, 1))) < 1e-10));

        vec o = trans(transition) * equilibrium;
        vec z = zeros(state_space);
        REQUIRE(approximately_equal(o, z, SMALL));

    }

    SECTION("simulate_sequence_over_phylogeny") {
        uvec s = Palantir::sample_sequence(equilibrium, 100);
        vector<SiteSimulation> sims = Simulate::sequence_over_phylogeny(
                tree, transition, sampling, s);
        REQUIRE(sims.size() == s.n_elem);
        for (auto s : sims) {
            REQUIRE(s.substitutions.size() == tree.n_nodes);
        }
    }
}

