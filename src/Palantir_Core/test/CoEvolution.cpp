#include "test.hpp"

TEST_CASE("CoEvolution")
{
    vec nucleotide_equilibrium({0.25, 0.25, 0.25, 0.25});
    vec fitness_1(test_fitness_1);
    vec fitness_2(test_fitness_2);
    GeneticCode g("Standard nuclear");
    mat delta(20, 20, fill::randu);

    mat nucleotide_transition = HasegawaKishinoYano::transition(
            nucleotide_equilibrium);

    SECTION("CoEvolution model works")
    {
        vec equilibrium = Palantir::CoEvolution::equilibrium(
                1000, 1e-8, nucleotide_equilibrium, fitness_1, fitness_2, delta, g);
        REQUIRE(equilibrium.size() == g.size * g.size);
        REQUIRE(sum(equilibrium) - 1.0 < SMALL);

        mat transition = Palantir::CoEvolution::transition(
                1000, 1e-8, nucleotide_transition, fitness_1, fitness_2, delta, g);

        vec z = zeros(equilibrium.size());
        vec rowsums = sum(transition, 1);
        REQUIRE(approximately_equal(rowsums, z, SMALL));

        vec o = transition * equilibrium;

        REQUIRE(approximately_equal(o, z, 1e-5));
        //cout << sum(abs(o - zeros(1, pi.size()))) << endl;

    }

}
