# M9 -- per-site interval tables reach the widget as a flat table, or not at all.
#
# simulate_with_shared_time_heterogeneity() gives every site its own mode
# history, so sim$intervals is a LIST of interval tables. RenderIntervals.js
# filters ONE flat table by node, so before the fix plot.Simulation() handed it
# the list unchanged, the filter matched nothing, and the interval track came out
# silently blank -- no error, no warning.
#
# The fix is on the R side only (the js contract is unchanged): a plot restricted
# to a single site is given that site's table, anything wider gets no interval
# track plus a warning naming the remedy.
#
# EXPECTED against the OLD library: FAIL -- plot(sim, sites = 1) puts a NESTED
#   array (one array of intervals per site) in the payload rather than a flat
#   table, and plot(sim) on several sites emits no warning.
# EXPECTED against the FIXED library: PASS -- plot(sim, sites = 1) puts site 1's
#   flat interval table in the payload, and the multi-site plot warns and sends
#   no interval table.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "M9"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

payload <- function(widget) jsonlite::fromJSON(as.character(widget$x))

# the js side treats an empty object as "no intervals" (Util.js is_empty_object),
# and jsonlite serializes R's NULL inside a list as {}
no_intervals <- function(p) is.null(p$intervals) ||
    (is.list(p$intervals) && !is.data.frame(p$intervals) && length(p$intervals) == 0)

result <- try({
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }
    tree <- mk("(((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0);")

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(11)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    switcher <- GeneralTimeReversible(equilibrium = c(.5, .5),
                                      exchangeability = matrix(c(0, 1, 1, 0), 2, 2))

    seq <- sample_sequence(model = m1, length = 3)
    sim <- simulate_with_shared_time_heterogeneity(
        phylogeny = tree, switching_model = switcher,
        substitution_models = list(m1, m2), sequence = seq, start_mode = 0,
        rate = 1, switching_rate = 3,
        segment_length = 0.05, tolerance = 1)   # tolerance 1 keeps the rescaler out

    if (!is.list(sim$intervals) || is.data.frame(sim$intervals) ||
        length(sim$intervals) != 3) {
        fail(sprintf("expected 3 per-site interval tables, got %s of length %d",
                     paste(class(sim$intervals), collapse = "/"),
                     length(sim$intervals)))
    }

    # 1. a single-site plot must send that site's FLAT table
    single <- payload(plot(sim, sites = 1))
    if (!is.data.frame(single$intervals)) {
        fail(sprintf("plot(sim, sites = 1) sent %s, not a flat interval table",
                     paste(class(single$intervals), collapse = "/")))
    }
    if (!all(c("node", "state", "from", "to") %in% names(single$intervals))) {
        fail(sprintf("flat interval table has columns %s",
                     paste(names(single$intervals), collapse = ",")))
    }
    if (nrow(single$intervals) != nrow(sim$intervals[[2]])) {
        fail(sprintf("payload has %d interval rows, site 1 has %d",
                     nrow(single$intervals), nrow(sim$intervals[[2]])))
    }
    if (!isTRUE(all.equal(single$intervals$node, sim$intervals[[2]]$node)) ||
        !isTRUE(all.equal(single$intervals$state, sim$intervals[[2]]$state))) {
        fail("the payload's interval table is not site 1's")
    }

    # 2. a multi-site plot must warn and send nothing
    warned <- character(0)
    multi <- withCallingHandlers(
        payload(plot(sim)),
        warning = function(w) {
            warned <<- c(warned, conditionMessage(w))
            invokeRestart("muffleWarning")
        })

    if (length(warned) == 0) {
        fail("plotting all 3 sites at once did not warn about per-site intervals")
    }
    if (!any(grepl("per-site intervals", warned, fixed = TRUE))) {
        fail(sprintf("unexpected warning(s): %s", paste(warned, collapse = "; ")))
    }
    if (!no_intervals(multi)) {
        fail("a multi-site plot still sent an interval structure to the widget")
    }

    # 3. a simulation with ONE shared table must be untouched
    mode_tree <- mk("(((A:0,B:0):0,(C:1,D:1):1):0);", type = "mode")
    flat_sim <- simulate_over_interval_phylogeny(
        tree, mode_tree, list(m1, m2), seq, start_mode = 0)
    flat <- payload(plot(flat_sim))
    if (!is.data.frame(flat$intervals) || nrow(flat$intervals) == 0) {
        fail("a simulation with a single shared interval table lost it")
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
