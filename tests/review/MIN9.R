# MIN9 -- subsetting an Alignment gives back an Alignment.
#
# An Alignment is a character matrix carrying a class and a `type` attribute. The
# default matrix `[` drops both, so alignment[, 1:3] came back as a plain
# character matrix: as.fasta() failed with "no applicable method for 'as.fasta'
# applied to an object of class 'matrix'", plot() would not dispatch, and
# .normalize() lost the `type` it switches on. A single row or column was dropped
# to a bare vector on top of that.
#
# EXPECTED against the OLD library: FAIL -- as.fasta(alignment[, 1:3], file)
#   errors.
# EXPECTED against the FIXED library: PASS -- the subset keeps class, `type` and
#   the taxon rownames, one-column and one-row slices stay matrices, an explicit
#   drop = TRUE still drops, and as.fasta() on the subset writes the right fasta.
#
# NOTE: `[.Alignment` is not matched by the NAMESPACE's exportPattern (the name
# starts with "["), so it needs an explicit S3method("[", Alignment) line, which
# is recorded in tests_review/HANDOFF.md. Without it this test fails with
# "no applicable method" exactly as the old library did.

suppressMessages(library(PalantiR))

ID <- "MIN9"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

result <- try({
    f <- tempfile(fileext = ".newick")
    writeLines("(((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0);", f)
    tree <- Phylogeny(f)

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(5)
    model <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    sim <- simulate_over_phylogeny(tree, model, sample_sequence(model = model, length = 5),
                                   rate = 1)
    alignment <- sim$alignment

    subset <- alignment[, 1:3]

    if (!inherits(subset, "Alignment")) {
        fail(sprintf("alignment[, 1:3] has class %s; `[.Alignment` is missing or not registered (see tests_review/HANDOFF.md)",
                     paste(class(subset), collapse = "/")))
    }
    if (!identical(attr(subset, "type"), attr(alignment, "type"))) {
        fail(sprintf("the subset's type is %s, expected %s",
                     paste(attr(subset, "type")), paste(attr(alignment, "type"))))
    }
    if (!identical(rownames(subset), rownames(alignment))) {
        fail("the subset lost its taxon names")
    }
    if (!identical(dim(subset), c(4L, 3L))) {
        fail(sprintf("the subset is %s, expected 4x3", paste(dim(subset), collapse = "x")))
    }
    if (!all(unclass(subset) == unclass(alignment)[, 1:3])) {
        fail("the subset does not hold the first three columns")
    }

    # the reported failure mode
    out <- tempfile(fileext = ".fasta")
    as.fasta(subset, out)
    lines <- readLines(out)
    if (length(lines) != 8) {
        fail(sprintf("as.fasta() on the subset wrote %d lines, expected 8",
                     length(lines)))
    }
    if (!identical(lines[1], ">A")) {
        fail(sprintf("first fasta line is %s, expected >A", lines[1]))
    }
    if (!identical(lines[2], paste(unclass(alignment)[1, 1:3], collapse = ""))) {
        fail("the fasta sequence is not the subset's first row")
    }

    # single row / single column must not be dropped to a vector
    one_column <- alignment[, 1]
    if (!inherits(one_column, "Alignment") || !identical(dim(one_column), c(4L, 1L))) {
        fail(sprintf("alignment[, 1] came back as %s %s",
                     paste(class(one_column), collapse = "/"),
                     paste(dim(one_column), collapse = "x")))
    }
    one_row <- alignment[1, ]
    if (!inherits(one_row, "Alignment") || !identical(dim(one_row), c(1L, 5L))) {
        fail(sprintf("alignment[1, ] came back as %s %s",
                     paste(class(one_row), collapse = "/"),
                     paste(dim(one_row), collapse = "x")))
    }

    # an explicit drop = TRUE must still drop
    dropped <- alignment[, 1, drop = TRUE]
    if (!is.null(dim(dropped)) || inherits(dropped, "Alignment")) {
        fail("drop = TRUE did not drop")
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
