#include "Simulate.hpp"

deque<double> Palantir::Simulate::poisson(double time, double rate)
{
    deque<double> times;
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
    deque<double> times;
    deque<ullong> states_from;
    deque<ullong> states_to;

    uniform_real_distribution<double> rnd_unif(0, 1);
    ullong c = start; // current state
    ullong j;

    // FIX (2026-08-17): the loop increment divided the exponential's RATE
    // PARAMETER by `rate` instead of dividing the deviate, giving a mean waiting
    // time of rate/-Q(c,c) rather than 1/(rate * -Q(c,c)) -- wrong by a factor of
    // rate^2 for every waiting time after the first. Site rate multipliers were
    // therefore compressed (see PROJECT-RECORD 8E.2: r = 4 gave 2.62x, not 4x).
    double t = rnd_exp(-transition.at(c, c)) / rate;
    for (; t <= time; t += rnd_exp(-transition.at(c, c)) / rate) {
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
// --- Exact time-change rescaling (opt-in alternative to the segment scheme).
// See the accompanying technical note: because every rescaled segment applies a
// scalar multiple of one generator, the out-of-equilibrium branch is exactly
// the homogeneous chain run for an intrinsic duration tau*. With p0 the
// forecast entering the branch and g the per-state class outflux, the budget
// delivered by intrinsic time tau is F(tau) = p0' (int_0^tau e^{Qu} du) g,
// strictly increasing with F'(tau) = p0' e^{Qtau} g -> 1, so F(tau*) = t*rate
// has a unique root. F and F' come from one exponential of the augmented
// matrix [[Q, I], [0, 0]].

static vec rescale_class_outflux(const mat& Q, const string& scaling_type,
                                 const Palantir::GeneticCode& g)
{
    vec out(Q.n_rows, fill::zeros);
    for(const Palantir::Codon& i : g) {
        for(const Palantir::Codon& j : g) {
            if(Palantir::Codon::_distance(i, j) <= 1 && i != j) {
                bool syn = Palantir::Codon::_synonymous(i, j);
                if(scaling_type == "substitution" ||
                   (scaling_type == "synonymous" && syn) ||
                   (scaling_type == "non-synonymous" && !syn)) {
                    out[i.index] += Q.at(i.index, j.index);
                }
            }
        }
    }
    return out;
}

static double rescale_solve_tau(const mat& Q, const mat& Estep_aug,
                                double dtau, const vec& p0, const vec& g,
                                double budget, vec& exit_forecast,
                                vec& knot_tau, vec& knot_F)
{
    // March the row vector [p', p' * int e^{Qu} du] through the precomputed
    // augmented step matrix. Each step costs one vector-matrix product; F and
    // f are exact at every grid point. Stop when the budget is bracketed and
    // invert the final interval with a monotone Hermite cubic.
    const uword n = Q.n_rows;
    rowvec v(2 * n, fill::zeros);
    v.head(n) = p0.t();
    rowvec v_prev = v;

    std::vector<double> taus(1, 0.0), Fs(1, 0.0);
    double f_prev = std::max(dot(p0, g), 1e-300);
    double F_prev = 0.0;
    double f_cur = f_prev, F_cur = 0.0;
    const ullong kmax = (ullong)(std::max(budget, 1.0) / dtau) * 100 + 1000;
    ullong k = 0;
    while (F_cur < budget && k < kmax) {
        v_prev = v;
        v = v * Estep_aug;
        k++;
        f_prev = f_cur; F_prev = F_cur;
        f_cur = std::max(as_scalar(v.head(n) * g), 1e-300);
        F_cur = as_scalar(v.tail(n) * g);
        taus.push_back(k * dtau);
        Fs.push_back(F_cur);
    }
    if (F_cur < budget) {
        throw logic_error("Exact rescaler failed to reach the event budget; "
                          "the class flux may be degenerate for this model.");
    }

    // Hermite inversion on the final interval: F has values (F_prev, F_cur)
    // and derivatives (f_prev, f_cur) at the ends; both positive, F monotone.
    double a = 0.5;
    for (int it = 0; it < 50; it++) {
        double h00 = (1 + 2 * a) * (1 - a) * (1 - a);
        double h10 = a * (1 - a) * (1 - a);
        double h01 = a * a * (3 - 2 * a);
        double h11 = a * a * (a - 1);
        double Fh = h00 * F_prev + h10 * dtau * f_prev
                  + h01 * F_cur + h11 * dtau * f_cur;
        double dh = 6 * a * (a - 1) * (F_prev - F_cur)
                  + (3 * a * a - 4 * a + 1) * dtau * f_prev
                  + (3 * a * a - 2 * a) * dtau * f_cur;
        double step = (Fh - budget) / std::max(dh, 1e-300);
        a -= step;
        if (a < 0) a = 0;
        if (a > 1) a = 1;
        if (std::fabs(Fh - budget) < 1e-12 * std::max(1.0, budget)) break;
    }
    double tau = (k - 1) * dtau + a * dtau;

    // Exit forecast: v_prev holds the exact state at the last full grid
    // point; one small exponential covers the partial step.
    exit_forecast = trans(rowvec(v_prev.head(n)) *
                          expmat(Q * (tau - (k - 1) * dtau)));

    knot_tau = conv_to<vec>::from(taus);
    knot_F = conv_to<vec>::from(Fs);
    ullong m = knot_tau.n_elem - 1;
    knot_tau[m] = tau;                    // final knot lands exactly on tau*
    knot_F[m] = budget;
    return tau;
}

// Piecewise-linear inverse of the tabulated time change: intrinsic time to
// delivered budget (monotone in both).
static double rescale_budget_at(const vec& knot_tau, const vec& knot_F, double tau_e)
{
    const uword K = knot_tau.n_elem - 1;
    if(tau_e <= 0) return 0.0;
    if(tau_e >= knot_tau[K]) return knot_F[K];
    // uniform grid except possibly the final (partial) interval
    double dtau = (K >= 2) ? knot_tau[1] : knot_tau[K];
    uword j = std::min((uword)(tau_e / dtau), K - 1);
    if(j + 1 == K && K >= 2 && tau_e < knot_tau[K - 1]) j = K - 2;
    double w = (tau_e - knot_tau[j]) /
               std::max(knot_tau[j + 1] - knot_tau[j], 1e-300);
    return knot_F[j] * (1.0 - w) + knot_F[j + 1] * w;
}

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
        string scaling_type,
        string rescale_method)
{
    if (rescale_method != "segments" && rescale_method != "exact") {
        throw logic_error("rescale_method must be \"segments\" or \"exact\"");
    }
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

    // Exact-mode inverse time-change tables, one entry per pushed segment
    // (empty vec when the segment is not an exact-mode stretch): events on
    // such a segment are reported at branch position start + F(tau_e)/rate,
    // keeping output identical in semantics to the segment method, including
    // on branches that mix rescaled and plain intervals.
    vector<deque<vec>> tree_tc_tau(tree.n_nodes);
    vector<deque<vec>> tree_tc_F(tree.n_nodes);
    vector<deque<double>> tree_tc_start(tree.n_nodes);

    // one augmented step matrix per mode for the exact method's marching
    const double EXACT_DTAU = 0.01;
    std::map<ullong, mat> exact_step;

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

                    // GUARD (2026-08-17): the segment rescaler renormalises Q
                    // against the transient distribution using
                    // MutationSelection::scaling. Two preconditions have to hold
                    // for that to mean anything, and neither was checked.
                    // See PROJECT-RECORD 8H.

                    // (1) MutationSelection::scaling sizes its accumulator by the
                    // number of sense codons, so it is only applicable to
                    // single-codon models. A CoEvolution ("codon_pair") model has
                    // n_codons^2 states and previously died here with an opaque
                    // Armadillo dimension error -- but only when the rescaler
                    // happened to engage, which for a 3600-state model it does not
                    // do at the default tolerance (8H.3).
                    if (local_Q[mode].n_rows != g.size) {
                        throw logic_error(
                            "The transient segment rescaler supports single-codon "
                            "models only (MutationSelection). This model has " +
                            to_string(local_Q[mode].n_rows) + " states against " +
                            to_string(g.size) + " sense codons in the active "
                            "genetic code. Either simulate without mode changes, "
                            "or disable the rescaler by passing a tolerance of 1 "
                            "or more.");
                    }

                    if (rescale_method == "exact" && scaling_type != "none") {
                        // Exact time change: one homogeneous stretch of the
                        // mode's own matrix for intrinsic duration tau*, with
                        // the budget (finish-start)*rate delivered identically.
                        // Event times on this branch are reported in intrinsic
                        // time, a strictly monotone reparameterisation of
                        // branch position.
                        vec g_class = rescale_class_outflux(local_Q[mode], scaling_type, g);
                        if (exact_step.count(mode) == 0) {
                            const uword nn = local_Q[mode].n_rows;
                            mat aug(2 * nn, 2 * nn, fill::zeros);
                            aug.submat(0, 0, nn - 1, nn - 1) = local_Q[mode];
                            aug.submat(0, nn, nn - 1, 2 * nn - 1) = eye(nn, nn);
                            exact_step[mode] = expmat(aug * EXACT_DTAU);
                        }
                        vec exit_forecast, knot_tau, knot_F;
                        double tau = rescale_solve_tau(local_Q[mode],
                                                       exact_step[mode],
                                                       EXACT_DTAU, current_pi,
                                                       g_class,
                                                       (finish - start) * rate,
                                                       exit_forecast,
                                                       knot_tau, knot_F);
                        tree_tc_tau[n].push_back(knot_tau);
                        tree_tc_F[n].push_back(knot_F);
                        tree_tc_start[n].push_back(start);
                        tree_segment_start[n].push_back(start);
                        tree_segment_end[n].push_back(start + tau / rate);
                        tree_pi[n].push_back(exit_forecast);
                        tree_Q[n].push_back(local_Q[mode]);
                        tree_S[n].push_back(local_S[mode]);
                        current_pi = exit_forecast;
                        current_Q = local_Q[mode];
                        continue;
                    }

                    // NOTE: no constraint on scaling_type. The rescaler must
                    // normalise the SAME quantity the models were built with --
                    // which the preceding commit guarantees by passing
                    // scaling_type through -- but WHICH quantity that is, is the
                    // caller's modelling choice, set by the units of the guide
                    // tree's branch lengths. An earlier revision of this guard
                    // required "synonymous"; that was wrong, because it assumed
                    // the invariant should be the Ne-independent rate rather than
                    // whatever the branch lengths denominate. With amino-acid
                    // branch lengths the correct invariant is non-synonymous.

                    // Iterate over small branch segments -
                    IntervalHistory segments(finish - start, segment_length);
                    segments.fast_forward(start);

                    for (ullong s = 0; s < segments.size; s++) {
                        double s_length = segments.time_to[s] - segments.time_from[s];
                        double pi_rmsd = rmsd(current_pi, local_pi[mode]);

                        if (pi_rmsd > tolerance) {
                            // rescaling on segment
                            mat I = eye(size(current_Q));
                            // FIX (2026-08-19): the transient forecast must advance at
                            // the same speed as the simulation it rescales. over_time
                            // divides waiting times by `rate`, i.e. the realised process
                            // runs at Q*rate, but this forecast advanced at Q alone, so
                            // for rate > 1 the schedule lagged the sequence it was
                            // normalising (and led it for rate < 1). The delivered
                            // number of substitutions per unit branch length during the
                            // transient then depended on the site's rate multiplier:
                            // measured on the population-shift sweep, stems delivered
                            // 4.3% too few substitutions at Ne 7,500 and 13.9% too few
                            // at Ne 17,000 (rate multipliers mean 1.20). Advancing the
                            // forecast by rate * s_length restores the invariant that
                            // one unit of branch length is one expected substitution of
                            // the scaled class, for every site.
                            current_pi = trans(current_pi.t() * ((current_Q * (s_length * rate)) + I));

                            double rho = MutationSelection::scaling(
                                    current_pi, local_Q[mode], scaling_type, g);

                            current_Q = local_Q[mode] / rho;
                            mat current_S = Palantir::sampling(current_Q);

                            tree_segment_start[n].push_back(segments.time_from[s]);
                            tree_segment_end[n].push_back(segments.time_to[s]);
                            tree_tc_tau[n].push_back(vec());
                            tree_tc_F[n].push_back(vec());
                            tree_tc_start[n].push_back(0.0);
                            tree_pi[n].push_back(current_pi);
                            tree_Q[n].push_back(current_Q);
                            tree_S[n].push_back(current_S);
                        } else {
                            // done rescaling - rest of branch

                            tree_segment_start[n].push_back(segments.time_from[s]);
                            tree_segment_end[n].push_back(node.length);
                            tree_tc_tau[n].push_back(vec());
                            tree_tc_F[n].push_back(vec());
                            tree_tc_start[n].push_back(0.0);
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
                    tree_tc_tau[n].push_back(vec());
                    tree_tc_F[n].push_back(vec());
                    tree_tc_start[n].push_back(0.0);
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
                        if (tree_tc_tau[n][i].n_elem > 0) {
                            // exact mode: map intrinsic event times to branch
                            // position so reporting matches the segment method
                            for (ullong e = 0; e < s.size; e++) {
                                double tau_e = s.time[e] * rate;
                                s.time[e] = tree_tc_start[n][i]
                                    + rescale_budget_at(tree_tc_tau[n][i],
                                                        tree_tc_F[n][i], tau_e) / rate;
                            }
                        } else {
                            s.fast_forward(i_start);
                        }
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
    vector<deque<double>> switch_times(tree.n_nodes);
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
                    vector<ullong> states = Palantir::Simulate::steps(
                            switching_sampling, start_state, n_switches);
                    deque<ullong> state(n_switches + 1);
                    deque<double> time_from(n_switches + 1);
                    deque<double> time_to(n_switches + 1);
                    state[0] = start_state;
                    time_from[0] = 0;
                    time_to[0] = switch_times[n][0];
                    for(ullong i = 1; i < n_switches; i++) {
                        state[i] = states[i];
                        time_from[i] = switch_times[n][i-1];
                        time_to[i] = switch_times[n][i];
                    }
                    ullong last = n_switches - 1;
                    if (node.length != switch_times[n][last]) {
                        state[n_switches] = states[last];
                        time_from[n_switches] = switch_times[n][last];
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
