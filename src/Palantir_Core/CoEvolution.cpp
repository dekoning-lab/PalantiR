#include "CoEvolution.hpp"

vec Palantir::CoEvolution::equilibrium(
        ullong population_size,
        double mutation_rate,
        const vec& nucleotide_equilibrium,
        const vec& fitness_1,
        const vec& fitness_2,
        const mat& delta,
        const GeneticCode& g)
{
    ullong n_pairs = g.size * g.size;
    const ullong& N = population_size;
    double max = numeric_limits<double>::min();
    vec equilibrium(n_pairs);

    double scale = 0;
    
    // Log-scale normalizing constant
    for(const Codon i : g) {
        for(const Codon j : g) {
            ullong aa_i = i.amino_acid;
            ullong aa_j = j.amino_acid;
            
            double selection = fitness_1[aa_i] + fitness_2[aa_j];
            double population_selection = 2.0 * N * selection * delta.at(aa_i, aa_j);
            if(max < population_selection) {
                max = population_selection;
            }
        }
    }
    
    for(const Codon i : g) {
        for(const Codon j : g) {
            ullong aa_i = i.amino_acid;
            ullong aa_j = j.amino_acid;
            ullong pair_index = CodonPair::index(i, j, g);
            
            double selection = fitness_1[aa_i] + fitness_2[aa_j];
            double population_selection = 2.0 * N * selection * delta.at(aa_i, aa_j);
            equilibrium.at(pair_index) = population_selection;
            
            double s = exp(population_selection - max);
            for (int k = 0; k < 3; k++) {
                ullong n_i = i.nucleotides[k];
                ullong n_j = j.nucleotides[k];
                equilibrium.at(pair_index) += (
                        log(nucleotide_equilibrium[n_i]) +
                        log(nucleotide_equilibrium[n_j]));
                s *= nucleotide_equilibrium[n_i] * nucleotide_equilibrium[n_j];
            }
            scale += s;
        }
    }
    scale = max + log(scale);
   
    equilibrium -= scale;
    equilibrium = exp(equilibrium);
    
    return equilibrium;
}

mat Palantir::CoEvolution::transition(
        ullong population_size,
        double mutation_rate,
        const mat& nucleotide_transition_rates,
        const vec& fitness_1,
        const vec& fitness_2,
        const mat& delta,
        const GeneticCode& g)
{
    ullong n_pairs = g.size * g.size;
    const ullong& N = population_size;
    mat transition(n_pairs, n_pairs, fill::zeros);
    
    // (i,j) -> (k, l)
    for(const Codon i : g) {
        for(const Codon j : g) {
            ullong c_ij = CodonPair::index(i, j, g);
            for(const Codon k : g) {
                for(const Codon l : g) {
                    ullong c_kl = CodonPair::index(k, l, g);
                    ullong ik_d = Codon::_distance(i, k);    // i -> k
                    ullong jl_d = Codon::_distance(j, l);    // j -> l
                    
                    if((ik_d == 1 && jl_d == 1) || (ik_d == 1 && jl_d == 0) || (ik_d == 0 && jl_d == 1)) {
                        double mu_ik = 1.0;
                        double mu_jl = 1.0;

                        if (ik_d == 1) {
                            pair<ullong, ullong> ik_t = Codon::_substitution(i, k);
                            mu_ik = nucleotide_transition_rates.at(ik_t.first, ik_t.second) * mutation_rate;
                        }
                        if (jl_d == 1) {
                            pair<ullong, ullong> jl_t = Codon::_substitution(j, l);
                            mu_jl = nucleotide_transition_rates.at(jl_t.first, jl_t.second) * mutation_rate;
                        }
                        
                        double joint_fitness_ij = (fitness_1[i.amino_acid] + fitness_2[j.amino_acid]) 
                            * delta.at(i.amino_acid, j.amino_acid);
                        double joint_fitness_kl = (fitness_1[k.amino_acid] + fitness_2[l.amino_acid]) 
                            * delta.at(k.amino_acid, l.amino_acid);

                        transition.at(c_ij, c_kl) = (2 * N * mu_ik * mu_jl *
                                MutationSelection::fixation_probability(N, joint_fitness_kl - joint_fitness_ij));
                        
                    }
                }
            }
        }
    }
    
    vec row_sums = sum(transition, 1);
    transition.diag() = -row_sums;
   
    return transition;
}

double Palantir::CoEvolution::scaling(
        const vec& equilibrium,
        const mat& transition,
        const string scaling_type,
        const GeneticCode& g)
{

    if (scaling_type == "none") {
        return sum(equilibrium);
    }

    vec scale(g.size * g.size, fill::zeros);
    // (i,j) -> (k, l)
    for(const Codon i : g) {
        for(const Codon j : g) {
            ullong c_ij = CodonPair::index(i, j, g);
            for(const Codon k : g) {
                for(const Codon l : g) {
                    ullong ik_d = Codon::_distance(i, k);    // i -> k
                    ullong jl_d = Codon::_distance(j, l);    // j -> l

                    // first site
                    if(ik_d <= 1 && i != k) {
                        if(scaling_type == "substitution") {
                            scale[c_ij] += transition.at(i.index, k.index);
                        } else if(scaling_type == "synonymous" && Codon::_synonymous(i, k)) {
                            scale[c_ij] += transition.at(i.index, k.index);
                        } else if(scaling_type == "non-synonymous" && (!Codon::_synonymous(i, k))) {
                            scale[c_ij] += transition.at(i.index, k.index);
                        }
                    }

                    // second site
                    if(jl_d <= 1 && j != l) {
                        if(scaling_type == "substitution") {
                            scale[c_ij] += transition.at(j.index, l.index);
                        } else if(scaling_type == "synonymous" && Codon::_synonymous(i, k)) {
                            scale[c_ij] += transition.at(j.index, l.index);
                        } else if(scaling_type == "non-synonymous" && (!Codon::_synonymous(i, k))) {
                            scale[c_ij] += transition.at(j.index, l.index);
                        }
                    }

                }
            }
        }
    }

    return sum(equilibrium % scale);
}
