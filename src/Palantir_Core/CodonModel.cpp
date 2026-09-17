#include "CodonModel.hpp"

string Palantir::CodonModel::canonical_scaling_type(const string& scaling_type)
{
    if(scaling_type == "standard") {
        return "substitution";
    }
    if(scaling_type == "synonymous" || scaling_type == "dS" ||
       scaling_type == "ds") {
        return "dS";
    }
    if(scaling_type == "none" || scaling_type == "substitution" ||
       scaling_type == "synonymous-per-codon" ||
       scaling_type == "non-synonymous") {
        return scaling_type;
    }
    throw logic_error("Unknown scaling type '" + scaling_type + "'. Valid "
                      "values are \"standard\" (an alias of \"substitution\"), "
                      "\"substitution\", \"synonymous\" (the dS gauge), "
                      "\"dS\", \"synonymous-per-codon\", "
                      "\"non-synonymous\" and \"none\".");
}

vec Palantir::CodonModel::class_outflux(
        const mat& transition,
        const string& scaling_type,
        const GeneticCode& g)
{
    const string st = canonical_scaling_type(scaling_type);
    if(transition.n_rows != g.size || transition.n_cols != g.size) {
        throw logic_error("Codon scaling requires a square single-codon rate "
                          "matrix with one row and column per sense codon");
    }

    vec outflux(g.size, fill::zeros);
    if(st == "none") {
        return outflux;
    }
    if(st == "dS") {
        throw logic_error("dS scaling requires a neutral reference generator; "
                          "use neutral_dS_scaling rather than class_outflux");
    }

    for(const Codon& i : g) {
        for(const Codon& j : g) {
            if(i == j || Codon::_distance(i, j) != 1) {
                continue;
            }
            const bool synonymous = Codon::_synonymous(i, j);
            if(st == "substitution" ||
               (st == "synonymous-per-codon" && synonymous) ||
               (st == "non-synonymous" && !synonymous)) {
                outflux[i.index] += transition.at(i.index, j.index);
            }
        }
    }
    return outflux;
}

double Palantir::CodonModel::scaling(
        const vec& equilibrium,
        const mat& transition,
        const string& scaling_type,
        const GeneticCode& g)
{
    const string st = canonical_scaling_type(scaling_type);
    if(equilibrium.n_elem != g.size) {
        throw logic_error("Codon scaling requires one equilibrium frequency "
                          "per sense codon");
    }
    if(st == "none") {
        return sum(equilibrium);
    }
    return sum(equilibrium % class_outflux(transition, st, g));
}

double Palantir::CodonModel::neutral_dS_scaling(
        const vec& neutral_equilibrium,
        const mat& neutral_transition,
        const GeneticCode& g)
{
    const double neutral_total = scaling(
        neutral_equilibrium, neutral_transition, "substitution", g);
    return neutral_total / 3.0;
}

vec Palantir::CodonModel::neutral_dS_outflux(
        const mat& neutral_transition,
        const GeneticCode& g)
{
    return class_outflux(neutral_transition, "substitution", g) / 3.0;
}
