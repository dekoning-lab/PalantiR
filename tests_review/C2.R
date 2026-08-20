# C2 -- switching_poisson gave interval i the mode states[i] instead of the mode
# states[i-1] that was actually entered at switch i-1. The first sampled mode was
# therefore thrown away, and the last two intervals of every branch always shared
# a mode (the tail interval correctly carries states[n_switches-1], which the
# off-by-one also handed to the interval before it).
#
# EXPECTED OUTCOME
#   Before the fix: FAIL. With a two-state switcher whose sampling matrix forces
#     a change of mode at every switch, the interval sequence on a branch must
#     alternate; instead intervals 0 and 1 share a mode, and so do the last two.
#   After the fix:  PASS. Consecutive intervals on a branch always differ in mode.
#
# The engine RNG is not seedable, so several replicates are run; the test fails if
# any branch in any replicate has two consecutive intervals in the same mode, and
# refuses to pass vacuously if no branch ever carried more than one interval.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "C2"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }

    tree <- mk("(((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0);")

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    # forced alternation: from either state the only reachable state is the other
    switcher <- GeneralTimeReversible(equilibrium = c(.5, .5),
                                      exchangeability = matrix(c(0, 1, 1, 0), 2, 2))

    n_rep <- 10
    max_intervals <- 0
    violations <- character(0)

    for (rep in seq_len(n_rep)) {
        s <- sample_sequence(model = m1, length = 2)
        sim <- simulate_with_shared_time_heterogeneity(
            phylogeny = tree, switching_model = switcher,
            substitution_models = list(m1, m2), sequence = s, start_mode = 0,
            rate = 1, switching_rate = 3,
            segment_length = 0.05, tolerance = 1)  # tolerance 1 keeps the rescaler out

        for (site in seq_along(sim$intervals)) {
            iv <- sim$intervals[[site]]
            for (nd in unique(iv$node)) {
                g <- iv[iv$node == nd, , drop = FALSE]
                max_intervals <- max(max_intervals, nrow(g))
                if (nrow(g) > 1 && any(diff(g$state) == 0)) {
                    violations <- c(violations, sprintf(
                        "rep %d site %d node %d: modes %s do not alternate",
                        rep, site, nd, paste(g$state, collapse = ",")))
                }
            }
        }
    }

    if (max_intervals < 2) {
        fail(sprintf("inconclusive: no branch carried >1 interval (max %d)",
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
