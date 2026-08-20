# M8 -- plot.Phylogeny takes the guide tree positionally and forwards `...`.
#
# Before the fix the method was declared plot.Phylogeny(x, ..., mode_phylogeny =
# NULL): `mode_phylogeny` could only be matched by name, so the documented (and
# error-message-recommended) positional call plot(tree, mode_tree) put the guide
# tree into `...`, and `...` was collected but never passed on to PhyloPlot().
# The result was a plain tree with no interval track, and every widget option a
# caller supplied was silently discarded.
#
# EXPECTED against the OLD library: FAIL -- the positional call yields
#   plot_intervals = false and no intervals in the payload, and show_labels =
#   FALSE never reaches the payload (it stays at the default true).
# EXPECTED against the FIXED library: PASS -- plot(tree, mode_tree) puts
#   plot_intervals = true and a non-empty interval table in the widget payload,
#   plot(tree, show_labels = FALSE) puts show_labels = false in it, and an
#   explicit plot_intervals = FALSE still overrides the branch default.

suppressMessages(library(PalantiR))

ID <- "M8"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

payload <- function(widget) jsonlite::fromJSON(as.character(widget$x))

result <- try({
    tree <- Phylogeny(system.file("extdata", "primates.newick", package = "PalantiR"))
    mode_tree <- Phylogeny(system.file("extdata", "primates_switch.newick",
                                       package = "PalantiR"), type = "mode")

    # 1. positional guide tree
    positional <- payload(plot(tree, mode_tree))

    if (!isTRUE(positional$options$plot_intervals)) {
        fail(sprintf("plot(tree, mode_tree) gave plot_intervals = %s; the guide tree was dropped",
                     paste(positional$options$plot_intervals, collapse = ",")))
    }
    if (!is.data.frame(positional$intervals) || nrow(positional$intervals) == 0) {
        fail("plot(tree, mode_tree) sent no interval table to the widget")
    }
    if (!all(c("node", "state", "from", "to") %in% names(positional$intervals))) {
        fail(sprintf("interval table has columns %s",
                     paste(names(positional$intervals), collapse = ",")))
    }

    # the named call must keep working too
    named <- payload(plot(tree, mode_phylogeny = mode_tree))
    if (!isTRUE(named$options$plot_intervals)) {
        fail("plot(tree, mode_phylogeny = mode_tree) no longer plots intervals")
    }

    # 2. `...` reaches PhyloPlot
    labelled <- payload(plot(tree, show_labels = FALSE))
    if (!identical(labelled$options$show_labels, FALSE)) {
        fail(sprintf("show_labels = FALSE did not reach the payload (got %s)",
                     paste(labelled$options$show_labels, collapse = ",")))
    }

    sized <- payload(plot(tree, line_width = 9, circle_size = 7))
    if (!identical(as.numeric(sized$options$line_width), 9) ||
        !identical(as.numeric(sized$options$circle_size), 7)) {
        fail("line_width / circle_size did not reach the payload")
    }

    # 3. an explicit plot_intervals must still win over the branch default
    overridden <- payload(plot(tree, mode_tree, plot_intervals = FALSE))
    if (!identical(overridden$options$plot_intervals, FALSE)) {
        fail("an explicit plot_intervals = FALSE was overridden by the branch default")
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
