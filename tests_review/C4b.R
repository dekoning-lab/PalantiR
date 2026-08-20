# C4b -- in exact mode a zero-length interval asks the time-change solver to
# deliver a budget of zero. The marching loop then never ran, k stayed at 0, and
# the subsequent (k - 1) * dtau wrapped on the unsigned k to ~1.8e17: the branch
# was simulated for an effectively infinite intrinsic duration and never returned.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL (timeout -- the child never comes back).
#   After the fix:  PASS. The zero-length branch is skipped, the simulation
#     returns well inside the time limit, and the rest of the tree is simulated.
#
# The simulation runs in a forked child with a wall-clock limit, because a C++
# loop of that length cannot be interrupted from R.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))
suppressMessages(library(parallel))

ID <- "C4b"
TIME_LIMIT <- 60
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

run_guarded <- function(fun, seconds) {
    job <- parallel::mcparallel(try(fun(), silent = TRUE))
    res <- parallel::mccollect(job, wait = FALSE, timeout = seconds)
    if (is.null(res)) {
        tools::pskill(job$pid, tools::SIGKILL)
        suppressWarnings(parallel::mccollect(job))
        return(list(ok = FALSE, why = sprintf("did not return within %g s", seconds)))
    }
    v <- res[[1]]
    if (is.null(v)) return(list(ok = FALSE, why = "child delivered no result (crashed)"))
    if (inherits(v, "try-error")) {
        return(list(ok = FALSE, why = gsub("[\r\n]+", " ", as.character(v))))
    }
    list(ok = TRUE, value = v)
}

result <- try({
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }

    # node 2 = A has length 0.0 and is assigned a mode its parent does not use
    tree <- mk("((A:0.0,B:1.0):1.0);")
    modes <- mk("((A:1.0,B:0.0):0.0);", "mode")

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    t0 <- Sys.time()
    r <- run_guarded(function() {
        s <- sample_sequence(model = m1, length = 5)
        sim <- simulate_over_interval_phylogeny(
            phylogeny = tree, mode_phylogeny = modes,
            models = list(m1, m2), sequence = s, start_mode = 0,
            rate = 1, segment_length = 0.01, tolerance = 0.001,
            rescale_method = "exact")
        d <- sim$substitutions
        list(n_zero_branch = sum(d$node == 2), n_other = sum(d$node != 2))
    }, seconds = TIME_LIMIT)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!r$ok) fail(sprintf("exact mode on a zero-length branch: %s", r$why))
    if (r$value$n_zero_branch != 0) {
        fail(sprintf("%d substitutions on a zero-length branch",
                     r$value$n_zero_branch))
    }
    if (r$value$n_other == 0) fail("no substitutions anywhere on the tree")
    if (elapsed > TIME_LIMIT) fail(sprintf("took %.1f s", elapsed))
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
