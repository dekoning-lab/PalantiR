#include "test.hpp"

TEST_CASE("equilibrium_to_fitness")
{
    vec fitness(test_fitness_1);
    vec equilibrium(test_equilibrium_1);

    vec f1 = Palantir::equilibrium_to_fitness(equilibrium, 1000);
    REQUIRE(approximately_equal(f1, fitness, 1e-3));
}


TEST_CASE("simulate_poisson")
{
    vector<double> a = Simulate::poisson(1, 10);
    vector<double> b = Simulate::poisson(1, 1);
    REQUIRE(a.size() > b.size());
}

TEST_CASE("sampling_matrix") 
{
    mat Q("-.75 .25 .25 .25; .25 -.75 .25 .25; .25 .25 -.75 .25; .25 .25 .25 -.75");
    mat s = Palantir::sampling(Q);
    REQUIRE(all(s.col(3) == 1));
}


TEST_CASE("simulate_over_time") 
{
    vec nucleotide_equilibruim("0.25 0.25 0.25 0.25");
    vec fitness(test_fitness_1);
    GeneticCode g("Standard nuclear");

    Phylogeny tree(test_tree);
    mat nucleotide_transition = HasegawaKishinoYano::transition(nucleotide_equilibruim);

    vec codon_equilibrium = MutationSelection::equilibrium(1000, 1e-8, nucleotide_equilibruim, fitness, g);
    mat codon_transition = MutationSelection::transition(1000, 1e-8, nucleotide_transition, fitness, g);
    double rho = MutationSelection::scaling(codon_equilibrium, codon_transition, "synonymous", g);
    codon_transition /= rho;
    mat codon_sampling = Palantir::sampling(codon_transition);

    REQUIRE(sum(codon_equilibrium) - 1.0 < SMALL);

    SECTION("sample_sequence") {
        uvec s= Palantir::sample_sequence(codon_equilibrium, 10);
    }

    SECTION("simulate_over_time") {


        SubstitutionHistory s = Simulate::over_time(codon_transition, codon_sampling);
        cout << s << endl;

        REQUIRE(is_sorted(s.time.begin(), s.time.end()));

        for (ullong i = 1; i < s.size; i++) {
            REQUIRE(s.state_from[i] == s.state_to[i - 1]);
        }
    }

    SECTION("simulate_over_phylogeny") {
        uvec s = Palantir::sample_sequence(codon_equilibrium, 1);
        SiteSimulation sim = Simulate::over_phylogeny(tree, codon_transition, codon_sampling, s[0]);
        REQUIRE(sim.substitutions.size() == tree.n_nodes);
        //root has no substitutions
        REQUIRE(sim.substitutions[0].size == 0);

    }

    SECTION("simulate_sequence_over_phylogeny") {
        uvec s = Palantir::sample_sequence(codon_equilibrium, 1000);
        vector<SiteSimulation> sims = Simulate::sequence_over_phylogeny(
                tree, codon_transition, codon_sampling, s);
        REQUIRE(sims.size() == s.size());
        for (auto s : sims) {
            REQUIRE(s.substitutions.size() == tree.n_nodes);
        }
    }

    SECTION("simulate_sequence_over_intervals") {
        Phylogeny tree(test_tree);
        // tree with modes ({0, 1}) instead of branch lengths
        Phylogeny mode_tree(test_tree_mode);
        vector<IntervalHistory> tree_intervals = tree.to_intervals(mode_tree);
        for(ullong i = 0; i < tree_intervals.size(); i++) {
            cout << "Node: " << i << endl;
            cout << tree_intervals[i] << endl;
        }

        uvec s = Palantir::sample_sequence(codon_equilibrium, 2);

        vector<vec> equilibrium;
        vector<mat> transition;
        vector<mat> sampling;
        uvec population_sizes({100, 1000});

        for(ullong i = 0; i < population_sizes.n_elem; i++){
            vec pi =  MutationSelection::equilibrium(population_sizes[i], 1e-8, nucleotide_equilibruim, fitness, g);
            mat Q = MutationSelection::transition(population_sizes[i], 1e-8, nucleotide_transition, fitness, g);
            double rho = MutationSelection::scaling(pi, Q, "synonymous", g);
            Q /= rho;
            mat S = Palantir::sampling(Q);

            equilibrium.push_back(pi);
            transition.push_back(Q);
            sampling.push_back(S);
        }

        Simulate::sequence_over_intervals(tree, tree_intervals, equilibrium, transition, sampling, s, 0, g);
    }

    SECTION("switching_poisson")
    {
        vec flat_pi(".5 .5");
        mat exch(".5 .5; .5 .5");
        mat gtr = GeneralTimeReversible::transition(flat_pi, exch);
        mat switching_sampling = Palantir::sampling(gtr);

        vector<vector<IntervalHistory>> s = Simulate::switching_poisson(tree, switching_sampling, 10, 0, 1);
        REQUIRE(s.size() == 10);
    }

    //SECTION("simulate_with_nested_heterogeneity") {
        //SubstitutionModel ms1 = MutationSelection::Model(10, 1e-8, hky, fitness, g);
        //SubstitutionModel ms2 = MutationSelection::Model(100, 1e-8, hky, fitness, g);
        //SubstitutionModel ms3 = MutationSelection::Model(1000, 1e-8, hky, fitness, g);
        //SubstitutionModel ms4 = MutationSelection::Model(10000, 1e-8, hky, fitness, g);

        //vec flat_pi("0.25 0.25 0.25 0.25");
        //mat exchangeability("0 0.2 0.1 0.3; 0.2 0 0.1 0.2; 0.1 0.1 0 0.1; 0.3 0.2 0.1 0");

        //SubstitutionModel gtr = GeneralTimeReversible::Model(flat_pi, exchangeability);
        //uvec seq = Palantir::sample_sequence(ms1, 10);

        //vector<SubstitutionModel> models({ms1, ms2, ms3, ms4});

        //auto sh = Palantir::simulate_over_phylogeny_with_nested_heterogeneity(tree, gtr, models, seq);
    //}
}
