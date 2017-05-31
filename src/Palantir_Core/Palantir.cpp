#include "Palantir.hpp"

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

vec Palantir::equilibrium_to_fitness(const vec& equilibrium, ullong population_size, double epsilon)
{
    vec pi(equilibrium);
    ullong N = 2 * population_size;
    uvec small = find(pi < epsilon);
    pi.elem(small).fill(epsilon);
    pi /= sum(pi);
    return (log(pi) - log(max(pi)) + N) / N;
}
