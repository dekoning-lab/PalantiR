#include "MutationSelection.hpp"

double Palantir::MutationSelection::fixation_probability(ullong population_size, double selection)
{
    if (selection == 0) {
        return 1.0 / (2.0 * population_size);
    }
    // FIX (2026-08-20): guard both exponential asymptotes. For selection below
    // about -709.78, exp(-selection) overflows so numerator and denominator are
    // both -Inf and the ratio is NaN. The true value, ~e^{(2N-1)s}, is at or
    // below the double underflow threshold there for any N >= 1, so return 0
    // explicitly. For large positive selection the denominator is exactly 1 in
    // double precision and the probability is the 1 - e^{-s} asymptote; return
    // it directly so the behaviour is explicit for either sign. Both cut-offs
    // sit where the unguarded formula still returns exactly the same value
    // (0, because the denominator overflows first, and -expm1(-s)), so results
    // in the normal range are unchanged.
    if (selection <= -709.0) {
        return 0.0;
    }
    if (selection >= 709.0) {
        return -expm1(-selection);
    }
    return (-expm1(-selection)) / (-expm1(-2.0 * population_size * selection));
}


vec Palantir::MutationSelection::equilibrium(
        ullong population_size,
        double mutation_rate,
        const vec& nucleotide_equilibrium,
        const vec& fitness,
        const GeneticCode& g)
{
    ullong n_codons = g.size;
    const ullong& N = population_size;
    // FIX (2026-08-20): the log-sum-exp offset was initialised to
    // numeric_limits<double>::min(), which is the smallest POSITIVE double
    // (2.2e-308), not the most negative one. With an all-negative fitness
    // vector no population_selection ever exceeded it, the offset stayed at
    // ~0 while every term underflowed to exp(-huge) = 0, and the equilibrium
    // came out all-NaN. Initialise to -infinity so the running maximum is
    // correct for any sign of fitness.
    double max = -numeric_limits<double>::infinity();
    vec codon_equilibrium(n_codons, fill::zeros);
    double scale = 0;

    // Log-scale normalizing constant
    for(const Codon& i : g) {
        double population_selection = 2.0 * N * fitness[i.amino_acid];
        if (max < population_selection) {
            max = population_selection;
        }
    }

    for(const Codon& i : g) {
        double population_selection = 2.0 * N * fitness[i.amino_acid];
        codon_equilibrium[i.index] = population_selection;
        double s = exp(population_selection - max);
        for(int k = 0; k < 3; k++) {
            ullong n_i = i.nucleotides[k];
            codon_equilibrium[i.index] += log(nucleotide_equilibrium[n_i]);
            s *= nucleotide_equilibrium[n_i];
        }
        scale += s;
    }
    scale = max + log(scale);

    codon_equilibrium -= scale;
    codon_equilibrium = exp(codon_equilibrium);
    // Should this be re-normalized?

    codon_equilibrium /= sum(codon_equilibrium);

    return codon_equilibrium;
}

mat Palantir::MutationSelection::transition(
        ullong population_size,
        double mutation_rate,
        const mat& nucleotide_transition,
        const vec& fitness,
        const GeneticCode& g)
{
    const ullong& N = population_size;
    mat codon_transition(g.size, g.size, fill::zeros);

    for(const Codon& i : g) {

        for(const Codon& j : g) {
            if(Codon::_distance(i, j) <= 1 && i != j) {

                pair<ullong, ullong> s = Codon::_substitution(i, j);
                double mu = nucleotide_transition.at(s.first,s.second) * mutation_rate;

                codon_transition.at(i.index, j.index) = (
                        2 * N * mu *
                        fixation_probability(N, fitness[j.amino_acid] - fitness[i.amino_acid]));
            }
        }
    }

    vec row_sums = sum(codon_transition, 1);
    codon_transition.diag() = -row_sums;

    return codon_transition;
}

double Palantir::MutationSelection::scaling(
        const vec& equilibrium,
        const mat& transition,
        const string scaling_type,
        const GeneticCode& g)
{

    string st = scaling_type;
    if (st != "none" && st != "substitution" && st != "synonymous" && st != "non-synonymous") {
        // Previously this printed to stdout and fell back to "none", i.e. no
        // scaling at all -- a silent wrong answer. The 2016 fork of this package
        // called the total-rate option "standard", so passing that name here
        // produced unscaled rate matrices with only a message on stdout, which
        // R does not show. Fail loudly instead.
        throw logic_error("Unknown scaling type '" + st + "'. Valid values are "
                          "\"none\", \"substitution\", \"synonymous\" and "
                          "\"non-synonymous\" (the 2016 fork's \"standard\" is "
                          "called \"substitution\" here).");
    }
    if (st == "none") {
        return sum(equilibrium);
    }
    vec scale(g.size, fill::zeros);

    for(const Codon& i : g) {
        for(const Codon& j : g) {
            if (Codon::_distance(i, j) <= 1 && i != j) {
                if(st == "substitution") {
                    scale[i.index] += transition.at(i.index, j.index);
                } else if(st == "synonymous" && Codon::_synonymous(i, j)) {
                    scale[i.index] += transition.at(i.index, j.index);
                } else if(st == "non-synonymous" && (!Codon::_synonymous(i, j))) {
                    scale[i.index] += transition.at(i.index, j.index);
                }
            }
        }
    }
    return sum(equilibrium % scale);
}

