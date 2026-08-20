#include "Palantir.hpp"
#include "GeneticCode.hpp"

// FIX (2026-08-20, C7): the single definition of the engine RNG. It used to be
// a `static mt19937 rng(rd());` in Palantir.hpp, which gave each translation
// unit a private, independently seeded copy. Seeded from random_device so that
// unseeded runs stay non-deterministic; set_palantir_seed() (R) / seed_rng()
// (C++) makes a run reproducible. NOTE: this engine is independent of R's RNG,
// so set.seed() does not affect it.
mt19937 Palantir::rng(random_device{}());

void Palantir::seed_rng(unsigned int seed)
{
    Palantir::rng.seed(seed);
}

bool approximately_equal(const vec& a, const vec& b, double tolerance)
{
    if(a.n_elem != b.n_elem) {
        return false;
    }
    vec diff = abs(a - b);
    return all(diff <= tolerance);
}

double rmsd(const vec& a, const vec& b)
{
    vec diff = a - b;
    return sqrt(mean(diff % diff));
}

double rnd_exp(double rate)
{
    exponential_distribution<double> exp(rate);
    return exp(Palantir::rng);
}

uvec Palantir::sample_sequence(const vec& equilibrium, ullong length)
{

    vec sampling_distribution = cumsum(equilibrium);
    uvec sequence(length);
    uniform_real_distribution<double> rnd_unif(0, 1);

    for (ullong i = 0; i < length; i++) {
        double u = rnd_unif(Palantir::rng);

        // search
        for (ullong j = 0; j < equilibrium.n_elem; j++) {
            if (u <= sampling_distribution[j]) {
                sequence[i] = j;
                break;
            }
        }
    }
    return sequence;
}

mat Palantir::sampling(const mat& transition)
{
    mat sampling(transition);

    for (ullong i = 0; i < transition.n_rows; i++) {
        double s = 0;
        for (ullong j = 0; j < transition.n_cols; j++) {
            if (i != j) {
                s += transition.at(i, j);
            }
            sampling.at(i, j) = s;
        }
        for (ullong j = 0; j < transition.n_cols; j++) {
            sampling.at(i, j) /= s;
        }
    }
    return sampling;
}

vec Palantir::equilibrium_to_fitness(const vec& equilibrium, double population_size, double epsilon)
{
    // FIX (2026-08-20): population_size is a real number. The formula is
    // smooth in N, and callers legitimately pass non-integer effective sizes
    // (e.g. a branch-length-weighted mean); the previous ullong parameter
    // silently truncated them.
    vec pi(equilibrium);
    double N = 2.0 * population_size;
    uvec small = find(pi < epsilon);
    pi.elem(small).fill(epsilon);
    pi /= sum(pi);
    return (log(pi) - log(max(pi)) + N) / N;
}

// FIX (2026-08-20, M2): mutational weight of each amino acid, i.e. the total
// probability its sense codons carry under the mutation process alone,
// m_A = sum_{c: aa(c) = A} prod_k nucleotide_equilibrium[c_k].
// This is exactly the factor the mutation-selection equilibrium multiplies
// exp(2 N psi_A) by when the codon distribution is collapsed to amino acids,
// and it is what equilibrium_to_fitness() has to divide out to hit a requested
// amino-acid distribution.
vec Palantir::amino_acid_mutational_weight(const vec& nucleotide_equilibrium, const GeneticCode& g)
{
    vec weight(AminoAcid::size, fill::zeros);

    for(const Codon& c : g) {
        double w = 1;
        for(int k = 0; k < 3; k++) {
            w *= nucleotide_equilibrium[c.nucleotides[k]];
        }
        weight[c.amino_acid] += w;
    }
    return weight;
}

// FIX (2026-08-20, M2): degeneracy-aware inversion, see the note in
// Palantir.hpp. psi_A = (log(pi_A) - log(m_A) - max(log(pi) - log(m)) + 2N)/2N
// gives exp(2 N psi_A) * m_A proportional to pi_A, so the amino-acid marginal
// of MutationSelection::equilibrium built on this fitness vector reproduces
// `equilibrium`. The offset keeps the convention of the two-argument form:
// the fittest amino acid has fitness exactly 1.
vec Palantir::equilibrium_to_fitness(
        const vec& equilibrium,
        double population_size,
        const vec& nucleotide_equilibrium,
        const GeneticCode& g,
        double epsilon)
{
    vec pi(equilibrium);
    double N = 2.0 * population_size;
    uvec small = find(pi < epsilon);
    pi.elem(small).fill(epsilon);
    pi /= sum(pi);

    vec weight = amino_acid_mutational_weight(nucleotide_equilibrium, g);
    if(weight.n_elem != pi.n_elem) {
        throw invalid_argument("equilibrium_to_fitness: `equilibrium` must have "
                               "one entry per amino acid");
    }
    // An amino acid with no sense codon under the active genetic code has
    // m_A = 0; clamping keeps the logarithm finite instead of producing +inf.
    uvec unreachable = find(weight < epsilon);
    weight.elem(unreachable).fill(epsilon);

    vec psi = log(pi) - log(weight);
    return (psi - max(psi) + N) / N;
}
