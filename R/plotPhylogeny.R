plot.Phylogeny <- function(phylogeny) PhyloPlot(phylogeny)

print.Phylogeny <- function(phylogeny) cat(phylogeny$string)

PhyloPlot <- function(
    phylogeny,
    substitutions = NULL,
    intervals = NULL,
    n_sites = NULL,
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
    line_width = 4) {

    tree <- fromJSON(phylogeny$json, simplifyVector = F)

    data <- toJSON(list(
        tree = tree,
        substitutions = substitutions,
        intervals = intervals,
        options = list(
            sim_type = sim_type,
            sites = sites,
            n_sites = n_sites,
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
            line_width = line_width)
    ), auto_unbox = T)

    htmlwidgets::createWidget(
        name = "PhyloPlot",
        x = data,
        width = width,
        height = height,
        package = "PalantiR",
        sizingPolicy = htmlwidgets::sizingPolicy(
            viewer.padding = 0,
            browser.fill = T
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
