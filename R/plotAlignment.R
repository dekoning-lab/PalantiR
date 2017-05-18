.aa_colors <- list("A"="#77dd88",
                  "G"="#77aa88",
                  "C"="#99ee66",
                  "D"="#55bb33",
                  "E"="#55bb33",
                  "N"="#55bb33",
                  "Q"="#55bb33",
                  "I"="#66bbff",
                  "L"="#66bbff",
                  "M"="#66bbff",
                  "V"="#66bbff",
                  "F"="#9999ff",
                  "W"="#9999ff",
                  "Y"="#9999ff",
                  "H"="#5555ff",
                  "K"="#ffcc77",
                  "R"="#ffcc77",
                  "P"="#eeaaaa",
                  "S"="#ff4455",
                  "T"="#ff4455")

AlignmentPlot <- function(alignment, width = NULL, height = NULL) {

    taxa <- row.names(alignment)
    colors <- apply(apply(alignment, 2, as_amino_acid), c(1,2), function(x) .aa_colors[[x]])

    sequences <- lapply(seq_len(nrow(alignment)), function(row) {
        list(index = row,
             taxon = taxa[[row]],
             sequence = data.frame(
                 state = alignment[row,],
                 color = colors[row, ]))
    })

    data <- jsonlite::toJSON(sequences, auto_unbox = T)

    htmlwidgets::createWidget(
        name = "AlignmentPlot",
        x = data,
        width = width,
        height = height,
        package = "PalantiR",
        sizingPolicy = htmlwidgets::sizingPolicy(
            browser.fill = T
        )
    )
}

AlignmentPlotOutput <- function(outputId, width = "100%", height = "100%"){
    htmlwidgets::shinyWidgetOutput(outputId, "AlignmentPlot", width, height, package = "PalantiR")
}

AlignmentPlotRender <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) expr <- substitute(expr) # force quoted
    htmlwidgets::shinyRenderWidget(expr, AlignmentPlotOutput, env, quoted = TRUE)
}

as.fasta.alignment <- function(alignemnt) {
    taxa <- row.names(alignment)
    for(row in seq_len(nrow(alignment))) {
        cat(">", taxa[row], "\n", sep = "")
        cat(alignment[row, ], "\n", sep = "")
    }
}
