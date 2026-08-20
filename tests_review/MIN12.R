# MIN12 -- click-to-download the whole svg is opt-in.
#
# PhyloPlot.js bound a click handler to the entire svg that downloaded it as a
# file whenever the widget was not in the RStudio viewer pane. Undocumented,
# unguarded and surprising in a browser or a Shiny app, where a click on the plot
# normally means something else. It is now behind a widget option,
# download_on_click, which PhyloPlot() emits and which defaults to FALSE.
#
# EXPECTED against the OLD library: FAIL -- there is no download_on_click option
#   in the payload and the js binds the handler unconditionally.
# EXPECTED against the FIXED library: PASS -- the option is in the payload,
#   false by default and true when asked for, and the js gates the handler on it.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "MIN12"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

options_of <- function(widget) jsonlite::fromJSON(as.character(widget$x))$options

result <- try({
    f <- tempfile(fileext = ".newick")
    writeLines("((A:1.0,B:1.0):1.0,C:1.0);", f)
    tree <- Phylogeny(f)

    # 1. the R side
    default <- options_of(plot(tree))
    if (is.na(match("download_on_click", names(default)))) {
        fail("PhyloPlot() does not emit a download_on_click option")
    }
    if (!identical(default$download_on_click, FALSE)) {
        fail(sprintf("download_on_click defaults to %s, expected FALSE",
                     paste(default$download_on_click, collapse = ",")))
    }

    asked <- options_of(PhyloPlot(tree, download_on_click = TRUE))
    if (!identical(asked$download_on_click, TRUE)) {
        fail("download_on_click = TRUE did not reach the payload")
    }

    # the option must also survive the plot() path, which forwards ...
    forwarded <- options_of(plot(tree, download_on_click = TRUE))
    if (!identical(forwarded$download_on_click, TRUE)) {
        fail("plot(tree, download_on_click = TRUE) did not reach the payload")
    }

    # 2. the js side must actually consult it (browser behaviour, so read the
    #    installed source; the marker is the guard around download_svg)
    path <- system.file("htmlwidgets", "PhyloPlot.js", package = "PalantiR")
    if (!nzchar(path)) {
        fail("PhyloPlot.js is not in the installed package")
    }
    js <- readLines(path, warn = FALSE)
    code <- js[!grepl("^\\s*//", js)]

    download <- grep("download_svg", code)
    if (length(download) == 0) {
        fail("PhyloPlot.js no longer contains the download handler at all")
    }
    guard <- grep("plot\\.options\\.download_on_click", code)
    if (length(guard) == 0) {
        fail("PhyloPlot.js binds the click handler without consulting download_on_click")
    }
    if (min(guard) > min(download)) {
        fail("the download_on_click guard comes after the download handler")
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
