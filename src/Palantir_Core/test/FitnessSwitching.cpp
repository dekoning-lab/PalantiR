#include "test.hpp"
TEST_CASE("FitnessSwitching")
{

    vec nucleotide_equilibrium("0.25 0.25 0.25 0.25");
    vec fitness_1(test_fitness_1);
    vec fitness_2(test_fitness_2);
    vector<vec> fitnesses({fitness_1, fitness_2});

    vec switching_equilibrium("0.5 0.5");

    GeneticCode g("Standard nuclear");

    Phylogeny tree(test_tree);
    mat exchangeability("0 0.1; 0.1 0");

    mat nucleotide_transition = HasegawaKishinoYano::transition(nucleotide_equilibrium);

    mat switching_transition = GeneralTimeReversible::transition(switching_equilibrium, exchangeability);

    ullong n_modes = fitnesses.size();
    ullong size = g.size;

    ullong state_space = size * n_modes;
    vector<vec> substitution_equilibrium;
    vector<mat> substitution_transition;
    vector<double> substitution_rho;

    for(ullong i = 0; i < n_modes; i++) {
        vec codon_equilibrium = MutationSelection::equilibrium(
                100, 1e-8, nucleotide_equilibrium, fitnesses[i], g);
        mat codon_transition = MutationSelection::transition(
                100, 1e-8, nucleotide_transition, fitnesses[i], g);
        double codon_rho = MutationSelection::scaling(
                codon_equilibrium, codon_transition, "synonymous", g);
        codon_transition /= codon_rho;

        substitution_equilibrium.push_back(codon_equilibrium);
        substitution_transition.push_back(codon_transition);
        substitution_rho.push_back(codon_rho);
    }

    mat transition = MarkovModulated::transition(
            substitution_transition, exchangeability, switching_equilibrium, g);

    vec equilibrium = MarkovModel::solve_equilibrium(transition);

    mat sampling = Palantir::sampling(transition);

    // end model construction

    SECTION("Switching Model check") {
        REQUIRE(sum(equilibrium) - 1.0 < SMALL);

        REQUIRE(all(vectorise(abs(sum(transition, 1))) < SMALL));

        vec o = trans(transition) * equilibrium;
        vec z = zeros(state_space);
        REQUIRE(approximately_equal(o, z, SMALL));

    }

    SECTION("simulate_sequence_over_phylogeny") {
        uvec s = Palantir::sample_sequence(equilibrium, 100);
        vector<SiteSimulation> sims = Simulate::sequence_over_phylogeny(
                tree, transition, sampling, s);

        REQUIRE(sims.size() == s.n_elem);
        for (auto s : sims) {
            REQUIRE(s.substitutions.size() == tree.n_nodes);
        }
    }
}

TEST_CASE("PoissonHeterogeneity")
{
    vec nucleotide_equilibrium("0.25 0.25 0.25 0.25");
    vec fitness_1(test_fitness_1);
    vec fitness_2(test_fitness_2);
    vector<vec> fitnesses({fitness_1, fitness_2});

    vec switching_equilibrium("0.5 0.5");

    GeneticCode g("Standard nuclear");

    Phylogeny tree(test_tree);
    mat exchangeability("0 0.1; 0.1 0");

    mat nucleotide_transition = HasegawaKishinoYano::transition(nucleotide_equilibrium);

    mat switching_transition = GeneralTimeReversible::transition(switching_equilibrium, exchangeability);
    mat switching_sampling = Palantir::sampling(switching_transition);


    vec pi_1 = MutationSelection::equilibrium(1000, 1e-8, nucleotide_equilibrium, fitness_1, g);
    vec pi_2 = MutationSelection::equilibrium(1000, 1e-8, nucleotide_equilibrium, fitness_2, g);
    mat Q_1 = MutationSelection::transition(1000, 1e-8, nucleotide_transition, fitness_1, g);
    mat Q_2 = MutationSelection::transition(1000, 1e-8, nucleotide_transition, fitness_2, g);
    double rho_1 = MutationSelection::scaling(pi_1, Q_1, "synonymous", g);
    double rho_2 = MutationSelection::scaling(pi_2, Q_2, "synonymous", g);
    Q_1 /= rho_1;
    Q_2 /= rho_2;
    mat S_1 = Palantir::sampling(Q_1);
    mat S_2 = Palantir::sampling(Q_2);

    vector<vec> equilibrium({pi_1, pi_2});
    vector<mat> transition({Q_1, Q_2});
    vector<mat> sampling({S_1, S_2});

    uvec sequence = Palantir::sample_sequence(pi_1, 10);
    ullong n_sites = sequence.n_elem;

    ullong n_modes = fitnesses.size();
    ullong start_mode = 0;
    ullong size = g.size;

    // Ludicrous speed
    double switching_rate = 100;
    double substitution_rate = 100;
    double segment_length = 0.0001;
    double tolerance = 0.001;

    vector<vector<Palantir::IntervalHistory>> tree_intervals = Palantir::Simulate::switching_poisson(
            tree, switching_sampling, n_sites, start_mode, switching_rate);

    vector<Palantir::SiteSimulation> sims;

    for(ullong site = 0; site < n_sites; site++) {
        cout << site << endl;
        uvec s(1);
        s[0] = sequence[site];

        vector<Palantir::SiteSimulation> sim = Palantir::Simulate::sequence_over_intervals(
             tree, tree_intervals[site], equilibrium, transition, sampling, s, start_mode, g, substitution_rate, segment_length, tolerance);

        REQUIRE(sim.size() == 1);
        sims.push_back(sim[0]);
    }
}
