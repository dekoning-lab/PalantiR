# C1 -- the segment rescaler's "done rescaling - rest of branch" segment used to
# end at node.length instead of at the finish of the interval it belongs to. On a
# branch carrying more than one mode interval that segment therefore ran to the
# end of the branch, and every later interval re-simulated on top of the overrun.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL. Substitution times on a multi-interval branch jump
#     backwards (a later interval starts before the previous one's overrun ends)
#     and the state chain recorded in production order is broken by those jumps.
#     Measured ~23% event inflation on the same construction.
#   After the fix:  PASS. Within every (site, node) group the recorded times are
#     non-decreasing, none exceeds the branch length, and every substitution
#     starts from the state the previous one ended in.
#
# The engine RNG is not seedable (set.seed has no effect on it), so the test runs
# several independent replicates and fails if ANY of them violates the invariant.
# It also refuses to pass vacuously: it checks that the construction really did
# produce multi-interval branches with rescaling engaged.

suppressMessages(library(PalantiR))

ID <- "C1"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }

    BRANCH <- 5.0
    tree <- mk(sprintf("(A:%.1f);", BRANCH))

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)                       # only fixes the fitness vectors
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    # two-state switcher: the sampling matrix forces a change of mode at every
    # switch, so a branch alternates modes and carries many intervals
    switcher <- GeneralTimeReversible(equilibrium = c(.5, .5),
                                      exchangeability = matrix(c(0, 1, 1, 0), 2, 2))

    n_rep <- 20
    max_intervals <- 0
    violations <- character(0)

    for (rep in seq_len(n_rep)) {
        s <- sample_sequence(model = m1, length = 2)
        sim <- simulate_with_shared_substitution_heterogeneity(
            phylogeny = tree, switching_model = switcher,
            substitution_models = list(m1, m2), sequence = s, start_mode = 0,
            rate = 1, switching_rate = 4,
            segment_length = 0.02, tolerance = 0.005)

        iv <- sim$intervals
        for (nd in unique(iv$node)) {
            max_intervals <- max(max_intervals, sum(iv$node == nd))
        }

        d <- sim$substitutions
        if (nrow(d) == 0) next
        for (st in unique(d$site)) for (nd in unique(d$node)) {
            g <- d[d$site == st & d$node == nd, , drop = FALSE]
            if (nrow(g) == 0) next
            if (any(diff(g$time) < 0)) {
                violations <- c(violations, sprintf(
                    "rep %d site %d node %d: event times not monotone", rep, st, nd))
            }
            if (max(g$time) > BRANCH + 1e-9) {
                violations <- c(violations, sprintf(
                    "rep %d site %d node %d: event at %.6g beyond branch length %.6g",
                    rep, st, nd, max(g$time), BRANCH))
            }
            if (nrow(g) > 1 &&
                !all(head(g$to, -1) == tail(g$from, -1))) {
                violations <- c(violations, sprintf(
                    "rep %d site %d node %d: state chain broken", rep, st, nd))
            }
        }
    }

    if (max_intervals < 3) {
        fail(sprintf("inconclusive: no branch carried >2 mode intervals (max %d)",
                     max_intervals))
    }
    if (length(violations) > 0) {
        fail(sprintf("%d violation(s), first: %s",
                     length(violations), violations[1]))
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
