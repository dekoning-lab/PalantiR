# C3 -- switching_poisson never wrote a branch's end-of-branch switching state
# back into tree_states, so every child branch restarted its mode chain from its
# parent's STARTING mode. The lineage mode process was not Markov across nodes.
#
# EXPECTED OUTCOME
#   Before the fix: FAIL. The first interval of every branch carries the root's
#     start_mode, whatever mode its parent branch ended in.
#   After the fix:  PASS. For every branch whose parent is itself a branch, the
#     first interval's mode equals the parent's last interval's mode.
#
# The engine RNG is not seedable, so several replicates are run. The test also
# refuses to pass vacuously: it requires that at least one parent branch actually
# ended in a mode other than the root's start_mode, which is the only situation in
# which the broken and fixed behaviours differ.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "C3"
START_MODE <- 0
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    # guide (mode) trees must be built with type = "mode"
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }

    tree <- mk("(((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0);")

    # parent map from the phylogeny's own json
    parent <- integer(0)
    walk <- function(node) {
        for (ch in node$children) {
            parent[[as.character(ch$index)]] <<- node$index
            walk(ch)
        }
    }
    walk(jsonlite::fromJSON(tree$json, simplifyVector = FALSE))
    root <- 0

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(7)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    switcher <- GeneralTimeReversible(equilibrium = c(.5, .5),
                                      exchangeability = matrix(c(0, 1, 1, 0), 2, 2))

    n_rep <- 10
    n_discriminating <- 0
    violations <- character(0)

    for (rep in seq_len(n_rep)) {
        s <- sample_sequence(model = m1, length = 2)
        sim <- simulate_with_shared_time_heterogeneity(
            phylogeny = tree, switching_model = switcher,
            substitution_models = list(m1, m2), sequence = s, start_mode = START_MODE,
            rate = 1, switching_rate = 3,
            segment_length = 0.05, tolerance = 1)

        for (site in seq_along(sim$intervals)) {
            iv <- sim$intervals[[site]]
            for (nd in unique(iv$node)) {
                p <- parent[[as.character(nd)]]
                if (is.null(p) || p == root) next   # parent is the root: no branch above
                child <- iv[iv$node == nd, , drop = FALSE]
                par <- iv[iv$node == p, , drop = FALSE]
                if (nrow(child) == 0 || nrow(par) == 0) next
                par_end <- par$state[nrow(par)]
                if (par_end != START_MODE) n_discriminating <- n_discriminating + 1
                if (child$state[1] != par_end) {
                    violations <- c(violations, sprintf(
                        "rep %d site %d node %d starts in mode %d but parent %d ends in mode %d",
                        rep, site, nd, child$state[1], p, par_end))
                }
            }
        }
    }

    if (n_discriminating == 0) {
        fail("inconclusive: no parent branch ended in a mode other than start_mode")
    }
    if (length(violations) > 0) {
        fail(sprintf("%d violation(s), first: %s",
                     length(violations), violations[1]))
    }
    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
