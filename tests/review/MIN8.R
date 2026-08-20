# MIN8 -- widget payloads are serialized at full precision.
#
# Every toJSON() call feeding a widget used jsonlite's default digits = 4, which
# rounds every number in the payload to four DECIMAL places. All tree geometry
# went through it, so a 2e-5 branch length was serialized as 0 and a shallow tree
# rendered collapsed onto its root.
#
# EXPECTED against the OLD library: FAIL -- the branch lengths in the payload of
#   a tree with 2e-5 branches are all 0.
# EXPECTED against the FIXED library: PASS -- they round-trip exactly, and a
#   length that four decimals cannot express (0.0001234) survives too.
#
# NOTE: the C++ newick-to-json writer itself emits 6 significant digits, so
# 0.000123456789 arrives in R already rounded to 0.000123457. That is a separate,
# upstream limit; this test stays inside it.

suppressMessages(library(PalantiR))

ID <- "MIN8"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

# the tree in the payload is a nested list of nodes; collect every `length`
collect_lengths <- function(node) {
    lengths <- if (is.null(node$length)) numeric(0) else as.numeric(node$length)
    for (child in node$children) {
        lengths <- c(lengths, collect_lengths(child))
    }
    lengths
}

result <- try({
    mk <- function(txt) {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f)
    }

    tiny <- mk("((A:0.00002,B:0.00002):0.00002,C:0.00002);")
    raw <- as.character(plot(tiny)$x)
    parsed <- jsonlite::fromJSON(raw, simplifyVector = FALSE)
    lengths <- collect_lengths(parsed$tree)   # the root's 0 is in here too

    if (length(lengths) == 0) {
        fail("the payload carries no branch lengths at all")
    }
    nonzero <- lengths[lengths > 0]
    if (length(nonzero) != 4) {
        fail(sprintf("expected 4 non-zero branch lengths in the payload, got %d (%s)",
                     length(nonzero), paste(lengths, collapse = ",")))
    }
    if (!isTRUE(all.equal(nonzero, rep(2e-05, 4)))) {
        fail(sprintf("branch lengths came back as %s, expected 2e-05",
                     paste(nonzero, collapse = ",")))
    }
    if (!grepl("2e-05", raw, fixed = TRUE)) {
        fail("the raw payload does not contain 2e-05")
    }

    # a length that four decimals cannot express at all (digits = 4 makes it 1e-04)
    fine <- mk("(A:0.0001234,B:0.0001234);")
    fine_lengths <- collect_lengths(
        jsonlite::fromJSON(as.character(plot(fine)$x), simplifyVector = FALSE)$tree)
    fine_lengths <- fine_lengths[fine_lengths > 0]
    if (length(fine_lengths) != 2 ||
        !isTRUE(all.equal(fine_lengths, rep(0.0001234, 2)))) {
        fail(sprintf("0.0001234 came back as %s",
                     paste(fine_lengths, collapse = ",")))
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
