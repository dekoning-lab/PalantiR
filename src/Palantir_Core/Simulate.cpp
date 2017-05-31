#include "Simulate.hpp"

vector<double> Palantir::Simulate::poisson(double time, double rate)
{
    vector<double> times;
    exponential_distribution<double> rnd_exp(rate);

    double t = rnd_exp(Palantir::rng);

    for (; t < time; t += rnd_exp(Palantir::rng)) {
        times.push_back(t);
    }
    return times;
}

vector<ullong> Palantir::Simulate::steps(
        const mat& sampling,
        ullong start,
        ullong steps)
{
    vector<ullong> sequence(steps);
    uniform_real_distribution<double> rnd_unif(0, 1);
    ullong c = start; // current state
    ullong j;

    for (ullong i = 0; i < steps; i++) {
        double r = rnd_unif(Palantir::rng);
        for (j = 0; r >= sampling(c, j); j++); // linear search
        sequence[i] = j;
        c = j;
    }
    return sequence;
}

Palantir::SubstitutionHistory Palantir::Simulate::over_time(
        const mat& transition,
        const mat& sampling,
        ullong start,
        double time,
        double rate)
{
    vector<double> times;
    vector<ullong> states_from;
    vector<ullong> states_to;

    uniform_real_distribution<double> rnd_unif(0, 1);
    ullong c = start; // current state
    ullong j;

    double t = rnd_exp(-transition.at(c, c)) / rate;
    for (; t <= time; t += rnd_exp(-transition.at(c, c) / rate)) {
        states_from.push_back(c); // starting state
        times.push_back(t);
        double r = rnd_unif(Palantir::rng);
        for (j = 0; r >= sampling.at(c, j); j++); // linear search
        states_to.push_back(j);
        c = j;
    }
    return SubstitutionHistory(times, states_from, states_to);
}

Palantir::SiteSimulation Palantir::Simulate::over_phylogeny(
        const Phylogeny& tree,
        const mat& transition,
        const mat& sampling,
        ullong start,
        double rate)
{
    uvec tree_states(tree.n_nodes);
    tree_states[0] = start;

    vector<reference_wrapper<const Phylogeny::Node> > nodes = tree.traversal();
    vector<SubstitutionHistory> tree_substitutions;

    for(const Phylogeny::Node& node : nodes) {
        if (!node.is_root()) {
            ullong n = node.index;
            ullong p = node.parent_index;
            tree_states[n] = tree_states[p];

            SubstitutionHistory s = Simulate::over_time(
                    transition,
                    sampling,
                    tree_states[n],
                    node.length,
                    rate);

            if (s.size) {
                tree_states[n] = s.state_to.back();
            }
            tree_substitutions.push_back(s);
        } else {
            tree_substitutions.push_back(SubstitutionHistory());
        }
    }
    return SiteSimulation(tree_substitutions, tree_states);
}

vector<Palantir::SiteSimulation> Palantir::Simulate::sequence_over_phylogeny(
        const Phylogeny& tree,
        const mat& transition,
        const mat& sampling,
        const uvec& sequence,
        double rate)
{
    vector<SiteSimulation> sims;

    for (ullong i = 0; i < sequence.n_elem; i++) {
        sims.push_back(Simulate::over_phylogeny(
                tree,
                transition,
                sampling,
                sequence[i],
                rate));
    }
    return sims;
}

// Shared modes for all sites in the sequence
vector<Palantir::SiteSimulation> Palantir::Simulate::sequence_over_intervals(
        const Phylogeny& tree,
        const vector<IntervalHistory>& tree_intervals,
        const vector<vec>& equilibrium,
        const vector<mat>& transition,
        const vector<mat>& sampling,
        const uvec& sequence,
        ullong start_mode,
        const GeneticCode& g,
        double rate,
        double segment_length,
        double tolerance,
        string scaling_type)
{
    if (tree_intervals.size() != tree.n_nodes) {
        throw logic_error("Intervals must correspond to tree nodes");
    }

    // Each node has multiple intervals
    vector<deque<vec>> tree_pi(tree.n_nodes);
    vector<vec> local_pi = equilibrium;
    tree_pi[0].push_back(local_pi[start_mode]);

    vector<deque<mat>> tree_Q(tree.n_nodes);
    vector<mat> local_Q = transition;
    tree_Q[0].push_back(local_Q[start_mode]);

    vector<deque<mat>> tree_S(tree.n_nodes);
    vector<mat> local_S = sampling;
    tree_S[0].push_back(local_S[start_mode]);

    vector<deque<double>> tree_segment_start(tree.n_nodes);

    vector<deque<double>> tree_segment_end(tree.n_nodes);

    vector<reference_wrapper<const Phylogeny::Node> > nodes = tree.traversal();

    // first find when rescaling is needed

    // for each branch (node)
    for(const Phylogeny::Node& node : nodes) {
        if (!node.is_root()) {
            ullong n = node.index;
            ullong p = node.parent_index;
            vec current_pi = tree_pi[p].back();
            mat current_Q = tree_Q[p].back();

            const IntervalHistory& intervals = tree_intervals[n];
            // for each interval on this node
            for (ullong i = 0; i < intervals.size; i++) {
                double start = intervals.time_from[i];
                double finish = intervals.time_to[i];
                ullong mode = intervals.state[i];

                double pi_rmsd = rmsd(current_pi, local_pi[mode]);

                if (pi_rmsd > tolerance) {
                    // rescaling

                    // Iterate over small branch segments -
                    IntervalHistory segments(finish - start, segment_length);
                    segments.fast_forward(start);

                    for (ullong s = 0; s < segments.size; s++) {
                        double s_length = segments.time_to[s] - segments.time_from[s];
                        double pi_rmsd = rmsd(current_pi, local_pi[mode]);

                        if (pi_rmsd > tolerance) {
                            // rescaling on segment
                            mat I = eye(size(current_Q));
                            current_pi = trans(current_pi.t() * ((current_Q * s_length) + I));

                            double rho = MutationSelection::scaling(
                                    current_pi, local_Q[mode], scaling_type, g);

                            current_Q = local_Q[mode] / rho;
                            mat current_S = Palantir::sampling(current_Q);

                            tree_segment_start[n].push_back(segments.time_from[s]);
                            tree_segment_end[n].push_back(segments.time_to[s]);
                            tree_pi[n].push_back(current_pi);
                            tree_Q[n].push_back(current_Q);
                            tree_S[n].push_back(current_S);
                        } else {
                            // done rescaling - rest of branch

                            tree_segment_start[n].push_back(segments.time_from[s]);
                            tree_segment_end[n].push_back(node.length);
                            tree_pi[n].push_back(local_pi[mode]);
                            tree_Q[n].push_back(local_Q[mode]);
                            tree_S[n].push_back(local_S[mode]);
                            break;
                        }
                    } // end rescaling segments
                } // end rescaling
                else { // entire interval
                    tree_segment_start[n].push_back(start);
                    tree_segment_end[n].push_back(finish);
                    tree_pi[n].push_back(local_pi[mode]);
                    tree_Q[n].push_back(local_Q[mode]);
                    tree_S[n].push_back(local_S[mode]);
                }
            } // end intervals
        } // end non root
    } // end nodes

    vector<SiteSimulation> sims;
    // simulate over sites
    for(ullong site = 0; site < sequence.n_elem; site ++) {

        vector<SubstitutionHistory> tree_substitutions;
        // Keep track of state and mode for each node
        uvec tree_states(tree.n_nodes);
        tree_states[0] = sequence[site];

        // simulate with tree_segments
        for(const Phylogeny::Node& node : nodes) {
            if (!node.is_root()) {
                ullong n = node.index;
                ullong p = node.parent_index;
                tree_states[n] = tree_states[p];
                SubstitutionHistory node_substitutions;

                for(ullong i = 0; i < tree_segment_start[n].size(); i++) {
                    double i_start = tree_segment_start[n][i];
                    double i_finish = tree_segment_end[n][i];

                    SubstitutionHistory s = Simulate::over_time(
                            tree_Q[n][i],
                            tree_S[n][i],
                            tree_states[n],
                            i_finish - i_start,
                            rate);

                    if (s.size) {
                        s.fast_forward(i_start);
                        tree_states[n] = s.state_to.back();
                        node_substitutions.append(s);
                    }
                } // end intervals
                tree_substitutions.push_back(node_substitutions);
            } // end root
            else {
                tree_substitutions.push_back(SubstitutionHistory());
            }
        } // end nodes

        sims.push_back(SiteSimulation(tree_substitutions, tree_states));
    } // end sites


    return sims;
}

vector<Palantir::IntervalHistory> Palantir::Simulate::switching_intervals(
        const Phylogeny& tree,
        const mat& switching_transition,
        const mat& switching_sampling,
        ullong start_mode,
        double rate)
{
    uvec tree_states(tree.n_nodes);
    tree_states[0] = start_mode;
    vector<reference_wrapper<const Phylogeny::Node> > nodes = tree.traversal();
    vector<IntervalHistory> tree_intervals(tree.n_nodes);
    for(const Phylogeny::Node& node : nodes) {
        if (!node.is_root()) {
            ullong n = node.index;
            ullong p = node.parent_index;
            tree_states[n] = tree_states[p];

            SubstitutionHistory s = Simulate::over_time(
                    switching_transition, switching_sampling,
                    tree_states[n], node.length, rate);
            if(s.size) {
                tree_intervals[n] = IntervalHistory(s, 0, node.length);
                tree_states[n] = s.state_to.back();
            } else {
                tree_intervals[n] = IntervalHistory(tree_states[n], 0, node.length);
            }
        }
    }
    return tree_intervals;
}

vector<vector<Palantir::IntervalHistory>> Palantir::Simulate::switching_poisson(
        const Phylogeny& tree,
        const mat& switching_sampling,
        ullong n_sites,
        ullong start_mode,
        double rate)
{
    vector<uvec> tree_states(n_sites);
    // Each site's root starts at the same state
    for(ullong site = 0; site < n_sites; site++) {
        tree_states[site] = uvec(tree.n_nodes);
        tree_states[site][0] = start_mode;
    }

    vector<reference_wrapper<const Phylogeny::Node> > nodes = tree.traversal();
    vector<vector<double>> switch_times(tree.n_nodes);
    vector<vector<IntervalHistory>> tree_intervals(n_sites);
    
    // Determine switch times for all sites
    for(const Phylogeny::Node& node : nodes) {
        if (!node.is_root()) {
            ullong n = node.index;
            switch_times[n] = Palantir::Simulate::poisson(node.length, rate);
        }
    }

    // For each site, simulate switches with times fixed
    for(ullong site = 0; site < n_sites; site++) {
        tree_intervals[site] = vector<IntervalHistory>(tree.n_nodes);
        for(const Phylogeny::Node& node : nodes) {
            if (!node.is_root()) {
                ullong n = node.index;
                ullong p = node.parent_index;
                tree_states[site][n] = tree_states[site][p];
                ullong n_switches = switch_times[n].size();
                if (n_switches) {
                    ullong start_state = tree_states[site][n];
                    vector<double> times = switch_times[n];
                    vector<ullong> states = Palantir::Simulate::steps(
                            switching_sampling, start_state, n_switches);
                    vector<ullong> state(n_switches + 1);
                    vector<double> time_from(n_switches + 1);
                    vector<double> time_to(n_switches + 1);
                    state[0] = start_state;
                    time_from[0] = 0;
                    time_to[0] = times[0];
                    for(ullong i = 1; i < n_switches; i++) {
                        state[i] = states[i];
                        time_from[i] = times[i-1];
                        time_to[i] = times[i];
                    }
                    ullong last = n_switches - 1;
                    if (node.length != times[last]) {
                        state[n_switches] = states[last];
                        time_from[n_switches] = times[last];
                        time_to[n_switches] = node.length;
                    }
                    tree_intervals[site][n] = IntervalHistory(state, time_from, time_to);
                } else {
                    tree_intervals[site][n] = IntervalHistory(tree_states[site][n], 0, node.length);
                }
            }
        }
    }

    return tree_intervals;
}
