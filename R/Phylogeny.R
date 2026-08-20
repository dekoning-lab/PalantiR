# FIX (2026-08-20, M8): two problems in one signature. `...` was collected and
# then never forwarded, so every widget option a caller passed to plot() --
# show_labels, sites, line_width, width/height -- was silently discarded. And
# `mode_phylogeny` was declared AFTER `...`, so it could only ever be matched by
# name: the documented positional call `plot(phylogeny, mode_phylogeny)` (which
# the error message below recommends) put the guide tree into `...`, where it
# was dropped, and quietly rendered a plain tree with no interval track.
plot.Phylogeny <- function(x, mode_phylogeny = NULL, ...) {
    phylogeny <- x
    if((phylogeny$type == "phylogeny") && is.null(mode_phylogeny)) {
       PhyloPlot(phylogeny, ...)
    } else if (phylogeny$type == "mode") {
        stop("To plot a mode phylogeny, use it as a second argument `plot(phylogeny, mode_phylogeny)`")
    } else if (
        (!is.null(mode_phylogeny)) &&
        (phylogeny$type == "phylogeny") &&
        (mode_phylogeny$type == "mode")) {
        intervals <- phylogeny_to_intervals(phylogeny, mode_phylogeny)
        # `plot_intervals` is the point of this branch, so it defaults to TRUE
        # here; do.call() rather than a direct call so that a caller who passes
        # it explicitly overrides the default instead of triggering "formal
        # argument matched by multiple actual arguments".
        options <- list(...)
        if(is.na(match("plot_intervals", names(options)))) {
            options$plot_intervals <- T
        }
        do.call(PhyloPlot, c(list(phylogeny, intervals = intervals), options))
    } else {
        stop("Unsupported phylogeny type")
    }
}

print.Phylogeny <- function(x, ...) cat(x$string)

PhyloPlot <- function(
    phylogeny,
    substitutions = NULL,
    intervals = NULL,
    width = NULL,
    height = NULL,
    sim_type = NULL,
    sites = "all",
    show_labels = T,
    plot_intervals = F,
    interval_opacity = 0.3,
    label_font_size = 10,
    show_axis = T,
    axis_ticks = 20,
    axis_tick_font_size = 10,
    padding = 10,
    circle_size = 4,
    circle_opacity = 1,
    line_width = 4,
    # FIX (2026-08-20, MIN12): clicking anywhere on the plot used to download the
    # whole svg whenever the widget was not in the RStudio viewer pane. That was
    # undocumented, unguarded and surprising in a browser or a Shiny app, so it
    # is now opt-in: set download_on_click = TRUE to get a click-to-save svg
    # (still suppressed inside the RStudio viewer, where the download does not
    # work).
    download_on_click = F) {

    tree <- fromJSON(phylogeny$json, simplifyVector = F)

    # FIX (2026-08-20, MIN8): toJSON()'s default digits = 4 rounds every number
    # in the payload to four decimal places, which silently flattened all tree
    # geometry -- a 2e-5 branch length was serialized as 0 and shallow trees
    # rendered collapsed onto the root. digits = NA serializes at full precision.
    data <- toJSON(list(
        tree = tree,
        substitutions = substitutions,
        intervals = intervals,
        options = list(
            sim_type = sim_type,
            sites = sites,
            show_labels = show_labels,
            plot_intervals = plot_intervals,
            interval_opacity = interval_opacity,
            label_font_size = label_font_size,
            show_axis = show_axis,
            axis_ticks = axis_ticks,
            axis_tick_font_size = axis_tick_font_size,
            padding = padding,
            circle_size = circle_size,
            circle_opacity = circle_opacity,
            line_width = line_width,
            download_on_click = download_on_click)
    ), auto_unbox = T, digits = NA)

    htmlwidgets::createWidget(
        name = "PhyloPlot",
        x = data,
        width = width,
        height = height,
        package = "PalantiR",
        sizingPolicy = htmlwidgets::sizingPolicy(
            viewer.padding = 0,
            viewer.fill = T,
            browser.fill = T,
            knitr.defaultHeight = 600,
            knitr.defaultWidth = 800
        )
    )
}

PhyloPlotOutput <- function(outputId, width = "100%", height = "400px"){
    htmlwidgets::shinyWidgetOutput(outputId, "PhyloPlot", width, height, package = "PalantiR")
}

PhyloPlotRender <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) expr <- substitute(expr) # force quoted
    htmlwidgets::shinyRenderWidget(expr, PhyloPlotOutput, env, quoted = TRUE)
}
