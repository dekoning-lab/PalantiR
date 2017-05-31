#include "test.hpp"

TEST_CASE("IntervalHistory") {
    vector<double> time({1, 2});
    vector<ullong> from({0, 1});
    vector<ullong> to({1, 2});
    SubstitutionHistory sh(time, from, to);
    IntervalHistory i(sh, 0, 3);

    REQUIRE(i.size == 3);

    REQUIRE(i.time_from[0] == 0);
    REQUIRE(i.time_to[0] == 1);
    REQUIRE(i.state[0] == 0);

    REQUIRE(i.time_from[1] == 1);
    REQUIRE(i.time_to[1] == 2);
    REQUIRE(i.state[1] == 1);

    REQUIRE(i.time_from[2] == 2);
    REQUIRE(i.time_to[2] == 3);
    REQUIRE(i.state[2] == 2);

    IntervalHistory j (SubstitutionHistory(), 0, 1);

    REQUIRE(j.size == 1);
    REQUIRE(j.state[0] == 0);
    REQUIRE(j.time_from[0] == 0);
    REQUIRE(j.time_to[0] == 1);

    IntervalHistory a(1.0, 0.3);
    REQUIRE(a.size == 4);
    REQUIRE(a.time_from[0] == 0);
    REQUIRE(a.time_to[0] == 0.3);
    REQUIRE(a.time_from[3] - 0.9 < TINY);
    REQUIRE(a.time_to[3] == 1.0);

    IntervalHistory b(2.0, 1.0);
    REQUIRE(b.size == 2);
    REQUIRE(b.time_from[0] == 0);
    REQUIRE(b.time_to[0] == 1);
    REQUIRE(b.time_from[1] == 1);
    REQUIRE(b.time_to[1] == 2);

    IntervalHistory c(1, 2);
    REQUIRE(c.size == 1);
    REQUIRE(c.time_from[0] == 0);
    REQUIRE(c.time_to[0] == 1);

}
