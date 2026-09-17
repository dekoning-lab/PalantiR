#include "test.hpp"
#include "../CodonModel.hpp"
#include "../GoldmanYang94.hpp"

TEST_CASE("GoldmanYang94")
{
    GeneticCode g("Standard nuclear");
    vec pi(g.size, fill::ones);
    pi /= sum(pi);
    const double omega = 0.25;
    const double kappa = 3.0;
    mat Q = GoldmanYang94::transition(pi, omega, kappa, g);

    SECTION("valid generator and stationary distribution")
    {
        REQUIRE(Q.n_rows == g.size);
        REQUIRE(Q.n_cols == g.size);
        REQUIRE(approximately_equal(sum(Q, 1), vec(g.size, fill::zeros), SMALL));
        REQUIRE(approximately_equal(trans(Q) * pi,
                                    vec(g.size, fill::zeros), SMALL));

        for(ullong i = 0; i < g.size; i++) {
            for(ullong j = 0; j < g.size; j++) {
                if(i != j) {
                    REQUIRE(std::abs(pi[i] * Q(i, j) - pi[j] * Q(j, i)) < SMALL);
                }
            }
        }
    }

    SECTION("omega and kappa multipliers")
    {
        ullong ttt = Codon::to_digit("TTT", g);
        ullong ttc = Codon::to_digit("TTC", g); // synonymous transition
        ullong tta = Codon::to_digit("TTA", g); // nonsynonymous transversion
        ullong ctt = Codon::to_digit("CTT", g); // nonsynonymous transition

        REQUIRE(std::abs(Q(ttt, ttc) / pi[ttc] - kappa) < SMALL);
        REQUIRE(std::abs(Q(ttt, tta) / pi[tta] - omega) < SMALL);
        REQUIRE(std::abs(Q(ttt, ctt) / pi[ctt] - kappa * omega) < SMALL);
    }

    SECTION("all codon scaling currencies")
    {
        const vector<string> types({"substitution", "synonymous-per-codon", "non-synonymous"});
        for(const string& type : types) {
            double rho = CodonModel::scaling(pi, Q, type, g);
            REQUIRE(rho > 0);
            REQUIRE(std::abs(CodonModel::scaling(pi, Q / rho, type, g) - 1.0) < SMALL);
        }
        REQUIRE(CodonModel::canonical_scaling_type("standard") == "substitution");
        REQUIRE(std::abs(CodonModel::scaling(pi, Q, "standard", g) -
                         CodonModel::scaling(pi, Q, "substitution", g)) < SMALL);
    }
}
