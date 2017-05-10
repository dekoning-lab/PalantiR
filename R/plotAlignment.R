AlignmentPlot <- function(alignment, width = NULL, height = NULL) {

    taxa <- dimnames(alignment)[[1]]

    sequences <- lapply(seq_len(nrow(alignment)), function(row) {
        list(index = row,
             sequence = alignment[row, ],
             taxon = taxa[[row]])
    })

    data <- jsonlite::toJSON(sequences)

    htmlwidgets::createWidget(
        name = "AlignmentPlot",
        x = data,
        width = width,
        height = height,
        package = "PalantiR"
    )
}

AlignmentPlotOutput <- function(outputId, width = "100%", height = "100%"){
    htmlwidgets::shinyWidgetOutput(outputId, "AlignmentPlot", width, height, package = "PalantiR")
}

AlignmentPlotRender <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) expr <- substitute(expr) # force quoted
    htmlwidgets::shinyRenderWidget(expr, AlignmentPlotOutput, env, quoted = TRUE)
}
