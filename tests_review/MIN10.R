# MIN10 -- AlignmentPlot clears its container before appending.
#
# AlignmentPlot.js's render() appended a fresh container div to the widget
# element without removing the previous one. htmlwidgets calls render() again on
# every Shiny re-render and on every resize, so each one stacked another complete
# copy of the alignment table below the last. RenderTree.js already clears its
# svg the same way (d3.select(...).selectAll("svg").remove()).
#
# This is browser-side behaviour that cannot be observed from R, so -- as for
# MIN2/MIN3/MIN4 -- the test asserts on the INSTALLED javascript: the clearing
# call must be present and must come before the container is appended. If the
# file is reformatted, update the markers grepped for here.
#
# EXPECTED against the OLD library: FAIL -- no remove() call in AlignmentPlot.js.
# EXPECTED against the FIXED library: PASS.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "MIN10"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    path <- system.file("htmlwidgets", "AlignmentPlot.js", package = "PalantiR")
    if (!nzchar(path)) {
        fail("AlignmentPlot.js is not in the installed package")
    }
    js <- readLines(path, warn = FALSE)

    clear <- grep("\\.remove\\(\\)", js)
    append_container <- grep("plot\\.container\\s*=\\s*d3\\.select", js)

    if (length(clear) == 0) {
        fail("AlignmentPlot.js never removes anything: a Shiny re-render stacks tables")
    }
    if (length(append_container) == 0) {
        fail("could not find the `plot.container = d3.select(...)` marker in AlignmentPlot.js")
    }
    if (!any(grepl("selectAll\\(\"div\"\\)|selectAll\\('div'\\)|selectAll\\(\"\\*\"\\)",
                   js[clear]))) {
        fail(sprintf("the remove() call does not clear the appended divs: %s",
                     trimws(js[clear[1]])))
    }
    if (min(clear) > min(append_container)) {
        fail("the container is appended before anything is removed")
    }

    # the widget itself must still build
    f <- tempfile(fileext = ".newick")
    writeLines("(((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0);", f)
    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(5)
    model <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    sim <- simulate_over_phylogeny(Phylogeny(f), model,
                                   sample_sequence(model = model, length = 4), rate = 1)
    widget <- plot(sim$alignment)
    if (!inherits(widget, "htmlwidget")) {
        fail("plot(alignment) no longer returns a widget")
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
