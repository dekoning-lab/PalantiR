#include "MutationSelection.hpp"

double Palantir::MutationSelection::fixation_probability(ullong population_size, double selection)
{
    if (selection == 0) {
        return 1.0 / (2.0 * population_size);
    } else{
        return (-expm1(-selection)) / (-expm1(-2.0 * population_size * selection));
    }
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
    double max = numeric_limits<double>::min();
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
        // need to raise error
        cout << "Unknown scaling type - defaulting to 'none'" << endl;
        st = "none";
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

