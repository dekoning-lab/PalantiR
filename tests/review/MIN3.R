# MIN3 -- rescale_budget_at carried a guard that re-seated the interpolation
# index j to K-2 when it landed on the final interval below knot_tau[K-1]. That
# combination cannot occur: landing on j == K-1 requires tau_e >= (K-1)*dtau ==
# knot_tau[K-1], which is the negation of the guard's own second condition. It is
# dead code and has been removed.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL (the dead guard is still in the source).
#   After the fix:  PASS.
#
# WHY THIS IS A SOURCE-LEVEL ASSERTION. Removing unreachable code is by
# definition invisible from outside; there is no input that takes the branch, so
# no behavioural test can distinguish the two versions. The test also exercises
# the exact rescaler end to end (which is what rescale_budget_at serves) to
# confirm the removal did not break the time-change inversion: event times must
# stay inside the branch and in order.

suppressMessages(library(PalantiR))

# Resolve the repository root by looking for src/Palantir_Core upward from the
# working directory; source-inspection tests SKIP when run without the source
# tree (for example against an installed package).
find_repo <- function() {
    for (cand in c("..", "../..", ".", "../../..")) {
        if (dir.exists(file.path(cand, "src", "Palantir_Core"))) return(normalizePath(cand))
    }
    NA
}

ID <- "MIN3"
repo <- find_repo()
if (is.na(repo)) { cat("MIN3: SKIP (source tree not available)\n"); quit(save = "no", status = 0) }
SRC <- file.path(repo, "src", "Palantir_Core", "Simulate.cpp")
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    if (!file.exists(SRC)) fail(sprintf("cannot read %s", SRC))
    src <- readLines(SRC, warn = FALSE)

    if (length(grep("rescale_budget_at", src, fixed = TRUE)) == 0) {
        fail("rescale_budget_at is gone from Simulate.cpp")
    }
    dead <- grep("j + 1 == K", src, fixed = TRUE)
    if (length(dead) > 0) {
        fail(sprintf("the dead j+1==K guard is still present at line %d", dead[1]))
    }

    # end-to-end check of the exact rescaler that rescale_budget_at serves
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }
    BRANCH <- 2.0
    tree  <- mk(sprintf("(A:%.1f);", BRANCH))
    modes <- mk("(A:1.0);", "mode")
    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    s <- sample_sequence(model = m1, length = 50)
    sim <- simulate_over_interval_phylogeny(
        tree, modes, list(m1, m2), s, 0,
        rate = 1, segment_length = 0.01, tolerance = 1e-9,
        rescale_method = "exact")
    d <- sim$substitutions
    if (nrow(d) == 0) fail("the exact rescaler produced no substitutions")
    if (any(!is.finite(d$time))) fail("the exact rescaler reported a non-finite event time")
    if (min(d$time) < 0 || max(d$time) > BRANCH + 1e-9) {
        fail(sprintf("the exact rescaler reported events outside [0, %.1f]: [%.6g, %.6g]",
                     BRANCH, min(d$time), max(d$time)))
    }
    for (st in unique(d$site)) {
        g <- d[d$site == st, , drop = FALSE]
        if (any(diff(g$time) < 0)) fail("exact-mode event times are not monotone")
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
