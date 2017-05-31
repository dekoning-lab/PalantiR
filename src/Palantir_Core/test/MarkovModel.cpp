#include "test.hpp"

TEST_CASE("MarkovModel")
{
    mat a = MarkovModel::random_transition(61);
    vec pi = MarkovModel::solve_equilibrium(a);

    vec o = trans(a) * pi;
    vec z = zeros(pi.n_elem);

    // ish
    REQUIRE(approximately_equal(o, z, 1e-5));
}
