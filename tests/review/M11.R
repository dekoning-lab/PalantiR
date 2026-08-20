# M11 -- the segmenting IntervalHistory constructor computed its size in closed
# form but filled the deques with an accumulating loop. On (end_time, length)
# pairs that are exact multiples in double arithmetic the accumulated position
# drifts just below the last grid point, so the loop ran one iteration too many:
# it wrote state[size] past the end of the deque, overwrote the final segment's
# start with the drifted end point, and thereby lost a whole segment of branch
# time. (0.8, 0.1) is such a pair: nothing at all was simulated in [0.7, 0.8).
#
# EXPECTED OUTCOME
#   Before the fix: FAIL. Over a branch of length 0.8 rescaled in segments of
#     0.1, the last tenth of the branch carries exactly zero substitutions --
#     deterministically, however many sites are simulated.
#   After the fix:  PASS. Every 0.1-wide window of the branch carries
#     substitutions and none is reported beyond the branch length.
#
# tolerance is set far below any attainable rmsd so that the rescaler runs every
# segment and never takes the "rest of branch" shortcut, which would otherwise
# paper over the lost segment.

suppressMessages(library(PalantiR))

ID <- "M11"
BRANCH <- 0.8
SEGMENT <- 0.1
N_SITES <- 200
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }

    tree  <- mk(sprintf("(A:%.1f);", BRANCH))
    modes <- mk("(A:1.0);", "mode")   # branch is mode 1, root starts in mode 0

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    s <- sample_sequence(model = m1, length = N_SITES)
    sim <- simulate_over_interval_phylogeny(
        phylogeny = tree, mode_phylogeny = modes,
        models = list(m1, m2), sequence = s, start_mode = 0,
        rate = 1, segment_length = SEGMENT, tolerance = 1e-9,
        rescale_method = "segments")

    d <- sim$substitutions
    if (nrow(d) == 0) fail("no substitutions simulated at all")
    if (max(d$time) > BRANCH + 1e-9) {
        fail(sprintf("substitution at %.6g beyond the branch length %.6g",
                     max(d$time), BRANCH))
    }

    edges <- seq(0, BRANCH, by = SEGMENT)
    counts <- as.integer(table(cut(d$time, breaks = edges, include.lowest = TRUE)))
    empty <- which(counts == 0)
    if (length(empty) > 0) {
        fail(sprintf("segment %d of %d ([%.2g, %.2g]) carries no substitutions out of %d; counts %s",
                     empty[1], length(counts), edges[empty[1]], edges[empty[1] + 1],
                     nrow(d), paste(counts, collapse = ",")))
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
