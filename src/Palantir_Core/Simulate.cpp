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

    // FIX (2026-08-20): a non-positive budget -- what a zero-length interval
    // asks for, and what a degenerate class flux can produce -- never entered
    // the marching loop, so k stayed at 0 and the `(k - 1) * dtau` below wrapped
    // on the unsigned k to ~1.8e17: an effectively infinite stretch of
    // simulation rather than an error or a no-op. Nothing is delivered over a
    // non-positive budget, so return an intrinsic duration of zero with the
    // entering forecast intact and a single-knot table, which rescale_budget_at
    // reads as "no budget anywhere".
    if (!(budget > 0)) {
        exit_forecast = p0;
        knot_tau = vec(1, fill::zeros);
        knot_F = vec(1, fill::zeros);
        return 0.0;
    }

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
    // Also the whole answer for a degenerate single-knot table (K == 0, from a
    // non-positive budget): knot_tau[0] is 0, so every positive tau_e leaves
    // here with a delivered budget of 0 and the K - 1 below is never reached.
    if(tau_e >= knot_tau[K]) return knot_F[K];
    // uniform grid except possibly the final (partial) interval
    double dtau = (K >= 2) ? knot_tau[1] : knot_tau[K];
    // FIX (2026-08-20): a guard here re-seated j to K - 2 when it landed on the
    // final interval below knot_tau[K - 1]. That cannot happen: reaching j ==
    // K - 1 requires tau_e >= (K - 1) * dtau == knot_tau[K - 1], which is the
    // negation of the guard's own second condition. Removed as dead code.
    uword j = std::min((uword)(tau_e / dtau), K - 1);
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
    // FIX (2026-08-20): "exact" used to fall through to the segment scheme
    // whenever the models carried no scaled class, so the caller's choice of
    // rescaler was silently ignored. The exact time change is defined by the
    // class flux it has to deliver; with no class there is nothing to solve
    // for. Say so rather than quietly running something else.
    if (rescale_method == "exact" && scaling_type == "none") {
        throw logic_error("rescale_method \"exact\" needs models built with a "
                          "scaled class: the exact time change solves for the "
                          "intrinsic duration that delivers the branch's budget "
                          "of scaled-class events, and scaling_type \"none\" "
                          "defines no such class. Rebuild the models with a "
                          "scaling_type, or use rescale_method \"segments\".");
    }
    if (tree_intervals.size() != tree.n_nodes) {
        throw logic_error("Intervals must correspond to tree nodes");
    }

    vector<vec> local_pi = equilibrium;
    vector<mat> local_Q = transition;
    vector<mat> local_S = sampling;

    // FIX (2026-08-20): the segment rescaler used to store a full Q and a full
    // sampling matrix for every segment it pushed -- ~60 KB per segment for a
    // 61-state model, measured at 3.3 GB peak on a single length-5 branch at the
    // shared simulators' default segment_length of 1e-4. Every rescaled segment
    // applies a scalar multiple of its own mode's generator, and Palantir::
    // sampling is invariant under that scalar, so a segment is fully described
    // by its mode and one number: over_time on (Q_mode, S_mode) at rate * scal
    // has exactly the exposure Q_mode * rate * scal * dt that the old
    // (Q_mode / rho, sampling(Q_mode / rho)) at rate had, with scal = 1 / rho.
    // Segments therefore carry a mode index and a scalar instead of matrices.
    vector<deque<ullong>> tree_segment_mode(tree.n_nodes);
    vector<deque<double>> tree_segment_scal(tree.n_nodes);
    vector<deque<double>> tree_segment_start(tree.n_nodes);
    vector<deque<double>> tree_segment_end(tree.n_nodes);

    // Exact-mode inverse time-change tables. One table per exact-mode stretch
    // lives in a shared pool and a segment records the index of its table, or
    // NO_TIME_CHANGE when it is an ordinary segment. Events on an exact stretch
    // are reported at branch position start + F(tau_e)/rate, keeping output
    // identical in semantics to the segment method, including on branches that
    // mix rescaled and plain intervals.
    const ullong NO_TIME_CHANGE = numeric_limits<ullong>::max();
    vector<vec> tc_tau;
    vector<vec> tc_F;
    vector<double> tc_start;
    vector<deque<ullong>> tree_segment_tc(tree.n_nodes);

    // What a branch hands to its children: the forecast where it ends and the
    // generator in force there, as the same (mode, scalar) pair. Only the last
    // entry of the old per-node deques of pi/Q/S was ever read, and a branch
    // that pushes no segment at all (see the zero-length interval guard below)
    // used to leave them empty for the child's .back() to dereference.
    vector<vec> tree_end_pi(tree.n_nodes);
    vector<ullong> tree_end_mode(tree.n_nodes, start_mode);
    vector<double> tree_end_scal(tree.n_nodes, 1.0);
    tree_end_pi[0] = local_pi[start_mode];

    // one augmented step matrix per mode for the exact method's marching
    const double EXACT_DTAU = 0.01;
    std::map<ullong, mat> exact_step;

    vector<reference_wrapper<const Phylogeny::Node> > nodes = tree.traversal();

    // Every segment push must write one entry to each of the parallel segment
    // tables. This is the only place that does so, so they cannot drift out of
    // alignment -- a misalignment silently reports events under the wrong
    // generator or at the wrong branch position.
    auto push_segment = [&](ullong index, double seg_start, double seg_end,
                            ullong seg_mode, double seg_scal, ullong seg_tc) {
        tree_segment_start[index].push_back(seg_start);
        tree_segment_end[index].push_back(seg_end);
        tree_segment_mode[index].push_back(seg_mode);
        tree_segment_scal[index].push_back(seg_scal);
        tree_segment_tc[index].push_back(seg_tc);
    };

    // first find when rescaling is needed

    // for each branch (node)
    for(const Phylogeny::Node& node : nodes) {
        if (!node.is_root()) {
            ullong n = node.index;
            ullong p = node.parent_index;
            vec current_pi = tree_end_pi[p];
            ullong current_mode = tree_end_mode[p];
            double current_scal = tree_end_scal[p];
            mat current_Q = local_Q[current_mode] * current_scal;

            const IntervalHistory& intervals = tree_intervals[n];
            // for each interval on this node
            for (ullong i = 0; i < intervals.size; i++) {
                double start = intervals.time_from[i];
                double finish = intervals.time_to[i];
                ullong mode = intervals.state[i];

                // FIX (2026-08-20): an interval of zero length -- every interval
                // on a zero-length branch, and any interval between coincident
                // switch times -- has nothing to simulate, but it used to be
                // handed to the segmenter anyway, which built a size-0
                // IntervalHistory and then wrote to its empty deques (segfault),
                // or in exact mode to the solver with a budget of zero (an
                // effectively infinite stretch, see rescale_solve_tau). Push no
                // segment; the child still inherits the forecast and generator
                // that entered the branch, because those are tracked below
                // rather than read off the last pushed segment.
                if (!(finish > start)) {
                    continue;
                }

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

                    if (rescale_method == "exact") {
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
                        tc_tau.push_back(knot_tau);
                        tc_F.push_back(knot_F);
                        tc_start.push_back(start);
                        push_segment(n, start, start + tau / rate,
                                     mode, 1.0, tc_tau.size() - 1);
                        current_pi = exit_forecast;
                        current_mode = mode;
                        current_scal = 1.0;
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

                            current_mode = mode;
                            current_scal = 1.0 / rho;
                            current_Q = local_Q[mode] * current_scal;

                            // the segmenter tiles [start, finish] exactly, but
                            // clamp anyway: no segment of an interval may cover
                            // branch time that belongs to the next one.
                            push_segment(n,
                                         std::min(segments.time_from[s], finish),
                                         std::min(segments.time_to[s], finish),
                                         current_mode, current_scal,
                                         NO_TIME_CHANGE);
                        } else {
                            // done rescaling - rest of branch

                            // FIX (2026-08-20): this segment used to end at
                            // node.length instead of at the finish of the
                            // interval it belongs to. On a branch carrying more
                            // than one mode interval it therefore ran past its
                            // own interval to the end of the branch, and every
                            // later interval re-simulated on top of the overrun:
                            // ~23% more events than the branch should carry
                            // (p = 3e-6), event times that jump backwards, and a
                            // state chain broken by the backwards jumps.
                            // FIX (2026-08-20): refresh the forecast and the
                            // generator here as well. They were left at the last
                            // rescaling segment's values, so a later rescaled
                            // interval on the same branch started its transient
                            // from a stale distribution.
                            current_pi = local_pi[mode];
                            current_mode = mode;
                            current_scal = 1.0;
                            current_Q = local_Q[mode];

                            push_segment(n,
                                         std::min(segments.time_from[s], finish),
                                         finish,
                                         current_mode, current_scal,
                                         NO_TIME_CHANGE);
                            break;
                        }
                    } // end rescaling segments
                } // end rescaling
                else { // entire interval
                    // FIX (2026-08-20): the forecast and generator were not
                    // refreshed on this path, so an interval that needed no
                    // rescaling left a later rescaled interval on the same
                    // branch starting its transient from whatever the previous
                    // interval happened to leave behind.
                    current_pi = local_pi[mode];
                    current_mode = mode;
                    current_scal = 1.0;
                    current_Q = local_Q[mode];

                    push_segment(n, start, finish, current_mode, current_scal,
                                 NO_TIME_CHANGE);
                }
            } // end intervals

            tree_end_pi[n] = current_pi;
            tree_end_mode[n] = current_mode;
            tree_end_scal[n] = current_scal;
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
                    ullong i_mode = tree_segment_mode[n][i];
                    ullong i_tc = tree_segment_tc[n][i];

                    // The segment's generator is its mode's, rescaled by one
                    // scalar; over_time divides its waiting times by the rate it
                    // is given, so folding the scalar into the rate delivers the
                    // same exposure Q * rate * scal * dt as a rescaled matrix
                    // would, and the sampling matrix does not depend on it. The
                    // scalar is exactly 1 on every non-rescaled segment.
                    SubstitutionHistory s = Simulate::over_time(
                            local_Q[i_mode],
                            local_S[i_mode],
                            tree_states[n],
                            i_finish - i_start,
                            rate * tree_segment_scal[n][i]);

                    if (s.size) {
                        if (i_tc != NO_TIME_CHANGE) {
                            // exact mode: map intrinsic event times to branch
                            // position so reporting matches the segment method
                            for (ullong e = 0; e < s.size; e++) {
                                double tau_e = s.time[e] * rate;
                                s.time[e] = tc_start[i_tc]
                                    + rescale_budget_at(tc_tau[i_tc],
                                                        tc_F[i_tc], tau_e) / rate;
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
                    // FIX (2026-08-20): interval i spans switch i-1 to switch i,
                    // so it carries the mode ENTERED at switch i-1, states[i-1].
                    // Indexing it with states[i] discarded the first sampled
                    // mode altogether and gave the last two intervals of every
                    // branch the same mode, since the tail interval below
                    // (correctly) also carries states[n_switches-1].
                    for(ullong i = 1; i < n_switches; i++) {
                        state[i] = states[i-1];
                        time_from[i] = switch_times[n][i-1];
                        time_to[i] = switch_times[n][i];
                    }
                    ullong last = n_switches - 1;
                    if (node.length != switch_times[n][last]) {
                        state[n_switches] = states[last];
                        time_from[n_switches] = switch_times[n][last];
                        time_to[n_switches] = node.length;
                    } else {
                        // no tail interval: drop the slot rather than leave a
                        // zero-width interval in mode 0 at the end of the branch
                        state.pop_back();
                        time_from.pop_back();
                        time_to.pop_back();
                    }
                    // FIX (2026-08-20): the mode the branch ends in was never
                    // written back, so every child restarted from its parent's
                    // STARTING mode and the lineage mode process was not Markov
                    // across nodes.
                    tree_states[site][n] = states[last];
                    tree_intervals[site][n] = IntervalHistory(state, time_from, time_to);
                } else {
                    tree_intervals[site][n] = IntervalHistory(tree_states[site][n], 0, node.length);
                }
            }
        }
    }

    return tree_intervals;
}
