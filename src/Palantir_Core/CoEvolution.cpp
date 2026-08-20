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
    // FIX (2026-08-20): the log-sum-exp offset was initialised to
    // numeric_limits<double>::min(), the smallest POSITIVE double, so an
    // all-negative fitness surface never updated it and the equilibrium came
    // out non-finite. Initialise to -infinity, as in MutationSelection.
    double max = -numeric_limits<double>::infinity();
    vec equilibrium(n_pairs);

    double scale = 0;
    
    // Log-scale normalizing constant
    for(const Codon& i : g) {
        for(const Codon& j : g) {
            ullong aa_i = i.amino_acid;
            ullong aa_j = j.amino_acid;
            
            double selection = fitness_1[aa_i] + fitness_2[aa_j];
            double population_selection = 2.0 * N * selection * delta.at(aa_i, aa_j);
            if(max < population_selection) {
                max = population_selection;
            }
        }
    }
    
    for(const Codon& i : g) {
        for(const Codon& j : g) {
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

    // FIX (2026-08-20): renormalise exactly as MutationSelection::equilibrium
    // does. The log-sum-exp constant leaves the vector summing to 1 only up to
    // floating-point error, and downstream code (scaling_type = "none",
    // simulation sampling) treats this as a probability distribution.
    equilibrium /= sum(equilibrium);

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
    for(const Codon& i : g) {
        for(const Codon& j : g) {
            ullong c_ij = CodonPair::index(i, j, g);
            for(const Codon& k : g) {
                for(const Codon& l : g) {
                    ullong c_kl = CodonPair::index(k, l, g);
                    ullong ik_d = Codon::_distance(i, k);    // i -> k
                    ullong jl_d = Codon::_distance(j, l);    // j -> l

                    // FIX (2026-08-20): simultaneous double substitutions
                    // (ik_d == 1 && jl_d == 1) previously got the rate
                    // 2N * mu_ik * mu_jl * fix -- the product of two mutation
                    // rates, which is dimensionally inconsistent (rate^2) and
                    // made the NORMALISED model depend on the absolute
                    // mutation_rate (54% of off-diagonal rate mass was double
                    // jumps at mutation_rate = 1, ~0% at 1e-6). In a standard
                    // CTMC two sites never change at the same instant; the
                    // instantaneous double rate is zero, and compensatory pairs
                    // arise as fast sequential single moves through the delta
                    // coupling.
                    if((ik_d == 1 && jl_d == 0) || (ik_d == 0 && jl_d == 1)) {
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
    if (scaling_type != "substitution" && scaling_type != "synonymous"
        && scaling_type != "non-synonymous") {
        // Without this an unrecognised name leaves the accumulator at zero and
        // the caller divides the transition matrix by it.
        throw logic_error("Unknown scaling type '" + scaling_type + "'. Valid "
                          "values are \"none\", \"substitution\", "
                          "\"synonymous\" and \"non-synonymous\".");
    }

    // FIX (2026-08-20): this function previously indexed the n_pairs x n_pairs
    // pair transition matrix with SINGLE-CODON indices (reading only the
    // top-left 61 x 61 corner, whose entries are pair-state rates, not
    // single-site rates), accumulated the first-site term inside the l loop and
    // the second-site term inside the k loop (a spurious x61 multiplicity), and
    // classified the second site with _synonymous(i, k) instead of (j, l). The
    // resulting time scale of every CoEvolution simulation was meaningless.
    //
    // Rewritten to walk the actual pair state space: for each pair state
    // P = (i, j) and target P' = (k, l) reachable at rate Q(P, P'), count how
    // many CODON changes of the move fall in the requested class -- a
    // single-site move contributes 1 if its changed codon is in the class, a
    // double move (rate now zero, see transition() above) would contribute 1
    // per changed codon in the class. The equilibrium-weighted total is the
    // expected number of class events per unit time per PAIR state; dividing
    // by 2 (two codon sites per pair state) makes one unit of branch length
    // equal one expected class event per CODON at equilibrium, exactly as in
    // MutationSelection::scaling. With delta == 1 (two independent sites) this
    // reduces to the mean of the two MutationSelection scaling factors.
    vec scale(g.size * g.size, fill::zeros);
    // (i,j) -> (k, l)
    for(const Codon& i : g) {
        for(const Codon& j : g) {
            ullong c_ij = CodonPair::index(i, j, g);
            for(const Codon& k : g) {
                for(const Codon& l : g) {
                    if(i == k && j == l) {
                        continue;
                    }
                    ullong c_kl = CodonPair::index(k, l, g);
                    double rate = transition.at(c_ij, c_kl);
                    if(rate == 0) {
                        continue;
                    }

                    double events = 0;
                    // first site: codon i -> k changed
                    if(i != k) {
                        if(scaling_type == "substitution"
                           || (scaling_type == "synonymous" && Codon::_synonymous(i, k))
                           || (scaling_type == "non-synonymous" && !Codon::_synonymous(i, k))) {
                            events += 1;
                        }
                    }
                    // second site: codon j -> l changed
                    if(j != l) {
                        if(scaling_type == "substitution"
                           || (scaling_type == "synonymous" && Codon::_synonymous(j, l))
                           || (scaling_type == "non-synonymous" && !Codon::_synonymous(j, l))) {
                            events += 1;
                        }
                    }

                    scale[c_ij] += rate * events;
                }
            }
        }
    }

    // Two codon sites per pair state: convert the per-pair-state event rate to
    // a per-codon-site event rate.
    return sum(equilibrium % scale) / 2.0;
}
