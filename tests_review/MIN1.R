# MIN1 -- rescale_method = "exact" was silently downgraded to the segment scheme
# whenever the models had been built with scaling_type = "none". The exact time
# change solves for the intrinsic duration that delivers a branch's budget of
# scaled-class events, so with no scaled class there is nothing for it to solve;
# quietly running a different rescaler instead hides that from the caller.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL. The call returns a simulation with no complaint.
#   After the fix:  PASS. The call raises an error naming scaling_type, and the
#     same models still simulate normally under rescale_method = "segments".

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "MIN1"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }

    tree  <- mk("(A:1.0);")
    modes <- mk("(A:1.0);", "mode")

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 0.5, hky, 1 + rnorm(20, 0, 5e-4),
                            scaling_type = "none")
    m2 <- MutationSelection(1000, 0.5, hky, 1 + rnorm(20, 0, 5e-4),
                            scaling_type = "none")
    s <- sample_sequence(model = m1, length = 5)

    thrown <- tryCatch({
        simulate_over_interval_phylogeny(
            tree, modes, list(m1, m2), s, 0,
            rate = 1, segment_length = 0.01, tolerance = 0.001,
            rescale_method = "exact")
        NULL
    }, error = function(e) conditionMessage(e))

    if (is.null(thrown)) {
        fail("rescale_method=\"exact\" with scaling_type=\"none\" returned a simulation instead of erroring")
    }
    if (!grepl("scaling_type", thrown, fixed = TRUE)) {
        fail(sprintf("error raised but it does not mention scaling_type: %s",
                     gsub("[\r\n]+", " ", thrown)))
    }

    # the same models must still work with the segment rescaler
    sim <- simulate_over_interval_phylogeny(
        tree, modes, list(m1, m2), s, 0,
        rate = 1, segment_length = 0.01, tolerance = 0.001,
        rescale_method = "segments")
    if (nrow(sim$substitutions) == 0) {
        fail("rescale_method=\"segments\" produced no substitutions for the same models")
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
