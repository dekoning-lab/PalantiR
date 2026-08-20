# M5 -- the guide tree must be a "mode" phylogeny and the primary tree a
# "phylogeny".
#
# Before the fix: the `type` field of a Phylogeny was never checked. Passing the
# branch-length tree in the `mode_tree` slot ran every branch under
# models[[1]] -- all modes collapsed to 0 -- with no warning, and passing the
# two trees the wrong way round was equally silent.
# EXPECTED against the OLD library: FAIL -- the mistyped trees are accepted.
# EXPECTED against the FIXED library: PASS -- both mistyped calls raise an error
#   naming the argument and the expected type, and the correctly typed call
#   still runs and uses more than one mode.
#
# NOTE: this makes `Phylogeny(guide, type = "mode")` mandatory for guide trees.
# Callers that relied on the default type must be updated (see
# tests_review/HANDOFF.md).

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

errs <- function(expr, must_mention) {
    m <- tryCatch({ force(expr); NA_character_ },
                  error = function(e) conditionMessage(e))
    if (is.na(m)) return(FALSE)
    all(vapply(must_mention, function(w) grepl(w, m, fixed = TRUE), logical(1)))
}

tryCatch({
    tree_path  <- system.file("extdata", "primates.newick", package = "PalantiR")
    guide_path <- system.file("extdata", "primates_switch.newick", package = "PalantiR")

    tree      <- Phylogeny(tree_path)
    guide     <- Phylogeny(guide_path, type = "mode")
    guide_bad <- Phylogeny(guide_path)              # default type "phylogeny"
    tree_mode <- Phylogeny(tree_path, type = "mode")

    hky <- HasegawaKishinoYano(rep(0.25, 4), transition_rate = 2)
    psi <- equilibrium_to_fitness(rep(1/20, 20), 1000)
    m   <- MutationSelection(1000, 1, hky, psi, scaling_type = "non-synonymous")
    models <- list(m, m, m)          # primates_switch.newick uses modes 1 and 2
    set_palantir_seed(7)
    seq <- sample_sequence(m, 2)

    sim_call <- function(a, b) simulate_over_interval_phylogeny(
        a, b, models, seq, start_mode = 1, rate = 1)

    checks <- c(
        "guide tree with type \"phylogeny\"" =
            errs(sim_call(tree, guide_bad), c("mode_tree", "mode")),
        "primary tree with type \"mode\"" =
            errs(sim_call(tree_mode, guide), c("tree", "phylogeny")),
        "phylogeny_to_intervals with a mistyped guide" =
            errs(phylogeny_to_intervals(tree, guide_bad), c("mode_phylogeny", "mode")),
        "Phylogeny() with an unknown type" =
            errs(Phylogeny(tree_path, type = "guide"), "type")
    )
    bad <- names(checks)[!checks]

    # Positive control: correctly typed trees still work and really do carry
    # more than one mode (the symptom the old silent path hid).
    iv <- phylogeny_to_intervals(tree, guide)
    sim <- sim_call(tree, guide)

    if (length(bad) > 0) {
        why <- paste("no error for:", paste(bad, collapse = "; "))
    } else if (length(unique(iv$state)) < 2) {
        why <- "positive control: the guide tree resolved to a single mode"
    } else if (is.null(sim$intervals) || nrow(sim$intervals) == 0) {
        why <- "positive control: simulation returned no intervals"
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "M5: PASS\n" else paste0("M5: FAIL ", why, "\n"))
