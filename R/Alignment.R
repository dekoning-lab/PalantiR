#http://www.bioinformatics.nl/~berndb/aacolour.html - MAEditor color scheme
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

#http://www.jalview.org/help/html/colourSchemes/
.nuc_colors <- list("A" = "#5CF659",
                    "C" = "#FFB14F",
                    "G" = "#F23C3F",
                    "T" = "#038BE9")

# Make sure the alignment is in "just codons"
#
# The result is a plain character matrix: everything below indexes it row by row
# and, since `[.Alignment` was added (MIN9), leaving the class on would hand
# those row slices back as one-row Alignments -- which toJSON() cannot serialize.
# The other two branches already return plain matrices from apply().
.normalize <- function(alignment, type) {

    if(type == "codon" || type == "nucleotide") {
        return(unclass(alignment));
    } else if(type == "compound_codon") {
        return(apply(alignment, c(1, 2), function(x)
            strsplit(x, ",")[[1]][[1]]))
    } else if(type == "codon_pair") {
        return(t(apply(alignment, 1, function(col)
            vapply(col, function(x)
                strsplit(x, ",")[[1]], character(2)))))
    } else {
        stop("Unknown alignment type")
    }
}

AlignmentPlot <- function(alignment, width = NULL, height = NULL) {
    type <- attr(alignment, "type")

    alignment <- .normalize(alignment, type)

    taxa <- row.names(alignment)

    if (grepl("codon", type)) {
        colors <- apply(apply(alignment, 2, as_amino_acid), c(1,2), function(x) .aa_colors[[x]])
    } else if (grepl("nucleotide", type)) {
        colors <- apply(alignment, c(1,2), function(x) .nuc_colors[[x]])
    }

    sequences <- lapply(seq_len(nrow(alignment)), function(row) {
        list(index = row,
             taxon = taxa[[row]],
             sequence = data.frame(
                 state = alignment[row,],
                 color = colors[row, ]))
    })

    # FIX (2026-08-20, MIN8): as in PhyloPlot(), toJSON()'s default digits = 4
    # rounds every number in the payload. Nothing here is numeric today, but the
    # two widget payloads should not disagree about precision.
    data <- jsonlite::toJSON(sequences, auto_unbox = T, digits = NA)

    htmlwidgets::createWidget(
        name = "AlignmentPlot",
        x = data,
        width = width,
        height = height,
        package = "PalantiR",
        sizingPolicy = htmlwidgets::sizingPolicy(
            viewer.padding = 0,
            browser.fill = T,
            viewer.fill = T,
            knitr.defaultHeight = 600,
            knitr.defaultWidth = 800
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

# FIX (2026-08-20, MIN9): an Alignment is a character matrix carrying a class and
# a `type` attribute, and the default matrix `[` drops both. Subsetting an
# alignment therefore produced a plain character matrix (or, for a single row or
# column, a bare vector), so alignment[, 1:3] could not be passed to as.fasta()
# ("no applicable method") or plot() and .normalize() lost the type it dispatches
# on. Restore class, `type` and the taxon rownames, and default to drop = FALSE
# so a one-site or one-taxon slice is still an alignment.
`[.Alignment` <- function(x, i, j, ..., drop = FALSE) {
    type <- attr(x, "type")
    subset <- NextMethod("[", drop = drop)
    if(is.matrix(subset)) {
        attr(subset, "type") <- type
        class(subset) <- "Alignment"
    }
    subset
}

as.fasta <- function(x, file = "", ...) UseMethod("as.fasta")
as.fasta.Alignment <- function(x, file = "", ...) {
    # plain matrix: `[.Alignment` (MIN9) would otherwise return the row slice
    # below as a one-row Alignment
    alignment <- unclass(x)
    taxa <- row.names(alignment)
    cat("", sep = "", file = file, append = F)
    for(row in seq_len(nrow(alignment))) {
        cat(">", taxa[row], "\n", sep = "", file = file, append = T)
        cat(alignment[row, ], "\n", sep = "", file = file, append = T)
    }
}

plot.Alignment <- function(x, ...) AlignmentPlot(x)
