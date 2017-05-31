#include "test.hpp"


TEST_CASE("MutationSelection")
{
    vec nucleotide_equilibrium("0.25 0.25 0.25 0.25");
    vec fitness(test_fitness_1);
    GeneticCode g("Standard nuclear");

    mat nucleotide_transition = HasegawaKishinoYano::transition(nucleotide_equilibrium);

    SECTION("Codon model check")
    {
        vec codon_equilibrium = MutationSelection::equilibrium(
                1000, 1e-8, nucleotide_equilibrium, fitness, g);

        mat codon_transition = MutationSelection::transition(
                1000, 1e-8, nucleotide_transition, fitness, g);

        double rho = MutationSelection::scaling(
                codon_equilibrium, codon_transition, "synonymous", g);

        codon_transition /= rho;

        REQUIRE(sum(codon_equilibrium) - 1.0 < SMALL);
        vec rowsums = sum(codon_transition, 1);
        vec z = zeros(g.size);
        REQUIRE(approximately_equal(rowsums, z, SMALL));

        vec o = trans(codon_transition) * codon_equilibrium;
        // ish
        REQUIRE(approximately_equal(o, z, 1e-4));
    }
}

