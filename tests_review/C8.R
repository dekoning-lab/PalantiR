# C8 -- fitness / delta dimensions are validated before construction.
#
# Before the fix: Models.cpp passed `fitness`, `fitness_1`, `fitness_2` and
# `delta` straight to kernels that index them with armadillo's UNCHECKED
# accessors. A 5-element fitness vector produced a MutationSelection model with
# 28 non-finite rates that differed between calls; a 2x2 delta produced a
# CoEvolution "model" whose equilibrium summed to 24.
# EXPECTED against the OLD library: FAIL -- the mis-sized arguments are accepted
#   and no error is raised.
# EXPECTED against the FIXED library: PASS -- each mis-sized argument raises an
#   error naming the argument, the expected size and the actual size, and a
#   correctly sized call still builds a finite model.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

# TRUE if expr raises an error whose message mentions all of `must_mention`.
errs <- function(expr, must_mention) {
    m <- tryCatch({ force(expr); NA_character_ },
                  error = function(e) conditionMessage(e))
    if (is.na(m)) return(FALSE)
    all(vapply(must_mention, function(w) grepl(w, m, fixed = TRUE), logical(1)))
}

tryCatch({
    hky  <- HasegawaKishinoYano(rep(0.25, 4), transition_rate = 2)
    good <- rep(1, 20)

    checks <- c(
        "MutationSelection fitness too short" =
            errs(MutationSelection(1000, 1, hky, rep(1, 5)), c("fitness", "20", "5")),
        "MutationSelection fitness too long" =
            errs(MutationSelection(1000, 1, hky, rep(1, 61)), c("fitness", "20", "61")),
        "CoEvolution fitness_1 wrong length" =
            errs(CoEvolution(100, 1, hky, rep(1, 19), good, matrix(1, 20, 20)),
                 c("fitness_1", "20", "19")),
        "CoEvolution fitness_2 wrong length" =
            errs(CoEvolution(100, 1, hky, good, rep(1, 21), matrix(1, 20, 20)),
                 c("fitness_2", "20", "21")),
        "CoEvolution delta 2x2" =
            errs(CoEvolution(100, 1, hky, good, good, matrix(1, 2, 2)),
                 c("delta", "20x20", "2x2")),
        "CoEvolution delta 20x4" =
            errs(CoEvolution(100, 1, hky, good, good, matrix(1, 20, 4)),
                 c("delta", "20x20", "20x4"))
    )

    bad <- names(checks)[!checks]

    # Positive control: a correctly sized fitness vector still builds a model
    # with an all-finite rate matrix.
    m <- MutationSelection(1000, 1, hky, equilibrium_to_fitness(rep(1/20, 20), 1000))
    if (length(bad) > 0) {
        why <- paste("no error for:", paste(bad, collapse = "; "))
    } else if (!all(is.finite(m$transition))) {
        why <- "positive control: valid 20-element fitness gave non-finite rates"
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "C8: PASS\n" else paste0("C8: FAIL ", why, "\n"))
