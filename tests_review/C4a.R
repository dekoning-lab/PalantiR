# C4a -- a zero-length branch whose mode differs from its parent's used to hand
# the segment rescaler an interval of length zero. The segmenting IntervalHistory
# constructor built a size-0 history and then wrote to its empty deques, which
# segfaults and takes the whole R session down.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL (the simulation crashes the child process, so nothing
#     comes back from it).
#   After the fix:  PASS. Both rescale methods complete on a tree with a
#     zero-length branch, the zero-length branch carries no substitutions, and the
#     rest of the tree is simulated normally.
#
# The simulation runs in a forked child so that a segfault or a hang is observed
# rather than inherited: a crash in the child leaves the parent free to report it.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))
suppressMessages(library(parallel))

ID <- "C4a"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

run_guarded <- function(fun, seconds) {
    job <- parallel::mcparallel(try(fun(), silent = TRUE))
    res <- parallel::mccollect(job, wait = FALSE, timeout = seconds)
    if (is.null(res)) {
        tools::pskill(job$pid, tools::SIGKILL)
        suppressWarnings(parallel::mccollect(job))
        return(list(ok = FALSE, why = "child died or timed out"))
    }
    v <- res[[1]]
    if (is.null(v)) {
        return(list(ok = FALSE, why = "child delivered no result (crashed)"))
    }
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

    # node 1 length 1.0 (mode 0), node 2 = A length 0.0 (mode 1), node 3 = B (mode 0)
    tree <- mk("((A:0.0,B:1.0):1.0);")
    modes <- mk("((A:1.0,B:0.0):0.0);", "mode")

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    for (method in c("segments", "exact")) {
        r <- run_guarded(function() {
            s <- sample_sequence(model = m1, length = 5)
            sim <- simulate_over_interval_phylogeny(
                phylogeny = tree, mode_phylogeny = modes,
                models = list(m1, m2), sequence = s, start_mode = 0,
                rate = 1, segment_length = 0.01, tolerance = 0.001,
                rescale_method = method)
            d <- sim$substitutions
            list(n = nrow(d),
                 n_zero_branch = sum(d$node == 2),
                 n_other = sum(d$node != 2))
        }, seconds = 120)

        if (!r$ok) fail(sprintf("rescale_method=%s: %s", method, r$why))
        if (r$value$n_zero_branch != 0) {
            fail(sprintf("rescale_method=%s: %d substitutions on a zero-length branch",
                         method, r$value$n_zero_branch))
        }
        if (r$value$n_other == 0) {
            fail(sprintf("rescale_method=%s: no substitutions anywhere on the tree",
                         method))
        }
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
