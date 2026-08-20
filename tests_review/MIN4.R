# MIN4 -- IntervalHistory.hpp declared the segmenting constructor as
# (length, end_time) while IntervalHistory.cpp defines it as (end_time, length).
# The names are the only documentation of the argument order, so anyone reading
# the header passed them the wrong way round.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL (the header's argument names are swapped relative to the
#     definition).
#   After the fix:  PASS.
#
# WHY THIS IS A SOURCE-LEVEL ASSERTION. Argument NAMES do not exist at runtime in
# C++, so the compiled library is identical either way and no behavioural test can
# distinguish them. The test also checks the constructor's live behaviour through
# the rescaler: the first argument really is the duration to cut up.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "MIN4"
HPP <- "/Users/jasondk/bot-workspace/incidental-convergence/palantir_pub/src/Palantir_Core/IntervalHistory.hpp"
CPP <- "/Users/jasondk/bot-workspace/incidental-convergence/palantir_pub/src/Palantir_Core/IntervalHistory.cpp"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

# Parameter names of a declaration/definition whose parameter list is exactly two
# `const double`s; NULL for any other signature.
double_pair <- function(lines, at) {
    txt <- ""
    for (l in lines[at:min(length(lines), at + 8)]) {
        txt <- paste(txt, l)
        if (grepl(")", l, fixed = TRUE)) break
    }
    if (!grepl("(", txt, fixed = TRUE) || !grepl(")", txt, fixed = TRUE)) return(NULL)
    inner <- sub("^[^(]*\\(", "", txt)
    inner <- sub("\\).*$", "", inner)
    args <- trimws(strsplit(inner, ",", fixed = TRUE)[[1]])
    if (length(args) != 2) return(NULL)
    if (!all(grepl("^const double [A-Za-z_][A-Za-z0-9_]*$", args))) return(NULL)
    sub("^const double ", "", args)
}

result <- try({
    if (!file.exists(HPP) || !file.exists(CPP)) fail("cannot read IntervalHistory sources")
    hpp <- readLines(HPP, warn = FALSE)
    cpp <- readLines(CPP, warn = FALSE)

    # the declaration whose two arguments are both plain doubles
    decl <- NULL
    for (at in grep("IntervalHistory\\($", hpp)) {
        a <- double_pair(hpp, at)
        if (length(a) == 2) { decl <- a; break }
    }
    if (is.null(decl)) fail("cannot find the two-double constructor declaration in the header")

    defn <- NULL
    for (at in grep("Palantir::IntervalHistory::IntervalHistory\\($", cpp)) {
        a <- double_pair(cpp, at)
        if (length(a) == 2) { defn <- a; break }
    }
    if (is.null(defn)) fail("cannot find the two-double constructor definition in the source")

    if (!identical(decl, defn)) {
        fail(sprintf("header declares (%s) but the definition reads (%s)",
                     paste(decl, collapse = ", "), paste(defn, collapse = ", ")))
    }
    if (!identical(defn, c("end_time", "length"))) {
        fail(sprintf("expected (end_time, length), found (%s)",
                     paste(defn, collapse = ", ")))
    }

    # behavioural cross-check: the rescaler cuts a branch of length 0.5 into
    # segments of 0.05, so events must span the whole branch and stop at it
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }
    BRANCH <- 0.5
    tree  <- mk(sprintf("(A:%.1f);", BRANCH))
    modes <- mk("(A:1.0);", "mode")
    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    s <- sample_sequence(model = m1, length = 200)
    sim <- simulate_over_interval_phylogeny(
        tree, modes, list(m1, m2), s, 0,
        rate = 1, segment_length = 0.05, tolerance = 1e-9,
        rescale_method = "segments")
    d <- sim$substitutions
    if (nrow(d) == 0) fail("no substitutions simulated")
    if (max(d$time) > BRANCH + 1e-9) {
        fail(sprintf("event at %.6g beyond the branch length %.1f", max(d$time), BRANCH))
    }
    if (max(d$time) < 0.9 * BRANCH) {
        fail(sprintf("events stop at %.6g, well short of the branch length %.1f",
                     max(d$time), BRANCH))
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
