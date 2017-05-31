#include "test.hpp"


TEST_CASE("HasegawaKishinoYano")
{
    vec nucleotide_equilibrium("0.25 0.25 0.25 0.25");
    mat nucleotide_transition = Palantir::HasegawaKishinoYano::transition(nucleotide_equilibrium);

    vec row_sums = sum(nucleotide_transition, 1);
    REQUIRE(sum(row_sums) == 0.0);
    REQUIRE(sum(nucleotide_equilibrium) == 1.0);

    vec o = nucleotide_transition * nucleotide_equilibrium;
    vec z = zeros(4);
    REQUIRE(approximately_equal(o, z, SMALL));
}
