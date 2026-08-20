# MIN11 -- the PhyloPlot tooltip no longer reads an option R never emits.
#
# Tooltip.js opened with
#     var site = (plot.options.n_sites === 1 ? 1 : d.site);
# `n_sites` is not in the options list PhyloPlot() builds, so the read was always
# undefined, and `site` was never used afterwards. Dead code that looked like a
# live per-site special case.
#
# Browser-side, so -- as for MIN2/MIN3/MIN4 -- the test asserts on the INSTALLED
# javascript plus a check that R still emits no such option.
#
# EXPECTED against the OLD library: FAIL -- n_sites appears in Tooltip.js.
# EXPECTED against the FIXED library: PASS -- it does not, no dead `site`
#   variable is left behind, and the tooltip's content function is intact.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "MIN11"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    path <- system.file("htmlwidgets", "PhyloPlot", "Tooltip.js", package = "PalantiR")
    if (!nzchar(path)) {
        fail("PhyloPlot/Tooltip.js is not in the installed package")
    }
    js <- readLines(path, warn = FALSE)
    code <- js[!grepl("^\\s*//", js)]          # ignore the comment explaining the fix

    if (any(grepl("n_sites", code))) {
        fail(sprintf("Tooltip.js still reads n_sites: %s",
                     trimws(code[grepl("n_sites", code)][1])))
    }
    if (any(grepl("var\\s+site\\s*=", code))) {
        fail("Tooltip.js still declares the unused `site` variable")
    }
    if (!any(grepl("var\\s+description\\s*=\\s*\\[name\\]", code))) {
        fail("Tooltip.js no longer builds its description from the tooltip name")
    }

    # and R must still not emit such an option
    f <- tempfile(fileext = ".newick")
    writeLines("((A:1.0,B:1.0):1.0,C:1.0);", f)
    options <- jsonlite::fromJSON(as.character(plot(Phylogeny(f))$x))$options
    if (!is.na(match("n_sites", names(options)))) {
        fail("PhyloPlot() emits an n_sites option after all")
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
