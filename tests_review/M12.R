# M12 -- the segment rescaler stored a full transition matrix AND a full sampling
# matrix for every segment it pushed (~60 KB per segment for a 61-state model;
# 3.3 GB peak on a single length-5 branch at the shared simulators' default
# segment_length of 1e-4). Every rescaled segment's generator is a scalar multiple
# of its mode's generator and the sampling matrix is invariant under that scalar,
# so a segment now carries a mode index and one scalar instead.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL on the memory check -- a length-3 branch rescaled in
#     segments of 1e-4 grows the process by well over a gigabyte (1.66 GB
#     measured here; the reviewer measured 3.3 GB at length 5).
#   After the fix:  PASS. The same branch costs a few megabytes, and the event
#     distribution is unchanged: a fully rescaled branch still delivers about one
#     scaled-class (synonymous) substitution per unit of branch length, under both
#     rescale methods.
#
# The distributional half of this test is a regression guard, not a
# before/after discriminator: it passes both before and after the refactor, and
# is here to catch a refactor that silently changed the sampled process.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))
suppressMessages(library(parallel))

ID <- "M12"
MEM_LIMIT_MB <- 300
N_SITES <- 800
BUDGET_LO <- 0.80   # ~5 s.e. at 800 sites; measured 0.9968 +- 0.0076 (segments)
BUDGET_HI <- 1.20   # and 0.9903 +- 0.0091 (exact) over 20 x 800-site replicates
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    # ---- 1. memory: a long branch cut into very short segments -------------
    big_tree  <- mk("(A:3.0);")
    big_modes <- mk("(A:1.0);", "mode")

    rss_kb <- function(pid) {
        x <- suppressWarnings(system2("ps", c("-o", "rss=", "-p", pid),
                                      stdout = TRUE, stderr = FALSE))
        if (length(x) == 0) NA_real_ else suppressWarnings(as.numeric(trimws(x[1])))
    }

    job <- parallel::mcparallel(try({
        s <- sample_sequence(model = m1, length = 1)
        sim <- simulate_over_interval_phylogeny(
            big_tree, big_modes, list(m1, m2), s, 0,
            rate = 1, segment_length = 1e-4, tolerance = 1e-9,
            rescale_method = "segments")
        nrow(sim$substitutions)
    }, silent = TRUE))

    base_kb <- NA_real_; peak_kb <- NA_real_; n_samples <- 0
    t0 <- Sys.time(); res <- NULL
    repeat {
        res <- parallel::mccollect(job, wait = FALSE, timeout = 0.02)
        v <- rss_kb(job$pid)
        if (!is.na(v)) {
            n_samples <- n_samples + 1
            if (is.na(base_kb)) base_kb <- v
            peak_kb <- if (is.na(peak_kb)) v else max(peak_kb, v)
        }
        if (!is.null(res)) break
        if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > 600) {
            tools::pskill(job$pid, tools::SIGKILL)
            suppressWarnings(parallel::mccollect(job))
            fail("the segment rescaler did not finish within 600 s")
        }
    }
    if (is.null(res[[1]]) || inherits(res[[1]], "try-error")) {
        fail(sprintf("the segment rescaler failed: %s",
                     if (is.null(res[[1]])) "child crashed"
                     else gsub("[\r\n]+", " ", as.character(res[[1]]))))
    }

    if (n_samples >= 2) {
        growth_mb <- (peak_kb - base_kb) / 1024
        if (growth_mb > MEM_LIMIT_MB) {
            fail(sprintf("a length-3 branch at segment_length 1e-4 grew the process by %.0f MB (limit %d MB)",
                         growth_mb, MEM_LIMIT_MB))
        }
    }

    # ---- 2. distribution: budget still delivered, both methods --------------
    tree  <- mk("(A:1.0);")
    modes <- mk("(A:1.0);", "mode")
    for (method in c("segments", "exact")) {
        s <- sample_sequence(model = m1, length = N_SITES)
        sim <- simulate_over_interval_phylogeny(
            tree, modes, list(m1, m2), s, 0,
            rate = 1, segment_length = 0.01, tolerance = 1e-9,
            rescale_method = method)
        d <- sim$substitutions
        per_site <- sum(d$synonymous) / N_SITES
        if (!is.finite(per_site) || per_site < BUDGET_LO || per_site > BUDGET_HI) {
            fail(sprintf("rescale_method=%s delivered %.3f synonymous substitutions per site over a branch of length 1 (expected ~1, band %.2f-%.2f)",
                         method, per_site, BUDGET_LO, BUDGET_HI))
        }
        if (max(d$time) > 1 + 1e-9) {
            fail(sprintf("rescale_method=%s reported an event at %.6g beyond the branch",
                         method, max(d$time)))
        }
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
