# MIN2 -- the transient forecast (current_pi) and the generator it is advanced
# with (current_Q) were not refreshed on the two paths that leave an interval
# without rescaling it: the "entire interval" path, and the "rest of branch"
# break that ends segment rescaling. A later rescaled interval on the same branch
# therefore started its transient from whatever the previous interval happened to
# leave behind.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL (the two paths do not refresh current_pi/current_Q).
#   After the fix:  PASS.
#
# WHY THIS IS A SOURCE-LEVEL ASSERTION. The numerical consequence is bounded by
# `tolerance`: both paths are taken precisely when the forecast is already within
# tolerance of the mode's equilibrium, so the stale and refreshed forecasts differ
# by at most that, and the delivered event counts differ by O(tolerance) times one
# segment. That is orders of magnitude below Monte-Carlo resolution on any
# affordable number of replicates, and the engine RNG cannot be seeded, so there
# is no honest behavioural discriminator to write. The test asserts the code
# change instead, and additionally checks that a branch mixing rescaled and
# non-rescaled intervals still produces coherent output.

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

ID <- "MIN2"
repo <- find_repo()
if (is.na(repo)) { cat("MIN2: SKIP (source tree not available)\n"); quit(save = "no", status = 0) }
SRC <- file.path(repo, "src", "Palantir_Core", "Simulate.cpp")
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    if (!file.exists(SRC)) fail(sprintf("cannot read %s", SRC))
    src <- readLines(SRC, warn = FALSE)

    block_after <- function(marker, n = 40) {
        at <- grep(marker, src, fixed = TRUE)
        if (length(at) != 1) return(NULL)
        src[at:min(length(src), at + n)]
    }

    refreshes <- function(block) {
        !is.null(block) &&
            any(grepl("current_pi = local_pi[mode];", block, fixed = TRUE)) &&
            any(grepl("current_Q = local_Q[mode];", block, fixed = TRUE))
    }

    rest <- block_after("// done rescaling - rest of branch")
    if (is.null(rest)) fail("cannot locate the \"rest of branch\" path in Simulate.cpp")
    if (!refreshes(rest)) {
        fail("the \"rest of branch\" path does not refresh current_pi/current_Q")
    }

    entire <- block_after("else { // entire interval")
    if (is.null(entire)) fail("cannot locate the \"entire interval\" path in Simulate.cpp")
    if (!refreshes(entire)) {
        fail("the \"entire interval\" path does not refresh current_pi/current_Q")
    }

    # behavioural sanity: a branch that mixes a non-rescaled interval with a
    # rescaled one still produces a coherent, in-bounds substitution record
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }
    BRANCH <- 4.0
    tree <- mk(sprintf("(A:%.1f);", BRANCH))
    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    switcher <- GeneralTimeReversible(equilibrium = c(.5, .5),
                                      exchangeability = matrix(c(0, 1, 1, 0), 2, 2))
    s <- sample_sequence(model = m1, length = 3)
    sim <- simulate_with_shared_substitution_heterogeneity(
        phylogeny = tree, switching_model = switcher,
        substitution_models = list(m1, m2), sequence = s, start_mode = 0,
        rate = 1, switching_rate = 3, segment_length = 0.02, tolerance = 0.005)
    d <- sim$substitutions
    if (nrow(d) == 0) fail("the mixed-interval branch produced no substitutions")
    for (st in unique(d$site)) {
        g <- d[d$site == st, , drop = FALSE]
        if (any(diff(g$time) < 0)) fail("event times on the mixed-interval branch are not monotone")
        if (max(g$time) > BRANCH + 1e-9) fail("an event was reported beyond the branch")
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
