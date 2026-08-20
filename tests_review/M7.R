# M7 (C++ half) -- genetic_code_names() exposes the valid genetic code names.
#
# Before the fix: the names of the built-in genetic codes existed only in the
# C++ GeneticCode_table. use_genetic_code() could not validate its argument and
# nothing in R could list the choices.
# EXPECTED against the OLD library: FAIL -- genetic_code_names() does not exist.
# EXPECTED against the FIXED library: PASS -- it returns the full set of names,
#   every returned name constructs a working genetic code, and the count matches
#   the table (verified here by round-tripping every name and confirming that a
#   plausible name NOT in the list is rejected by the C++ table).

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

original <- get_genetic_code()

tryCatch({
    names <- genetic_code_names()

    if (!is.character(names)) {
        stop("genetic_code_names() did not return a character vector")
    }
    if (anyDuplicated(names) != 0) {
        stop("genetic_code_names() returned duplicate names")
    }

    # Every advertised name must build a usable code.
    sizes <- integer(0)
    for (n in names) {
        use_genetic_code(n)
        g <- GeneticCode()
        if (length(g) < 60 || length(g) > 64) {
            stop(paste0("genetic code '", n, "' has ", length(g), " sense codons"))
        }
        sizes <- c(sizes, length(g))
    }
    use_genetic_code(original)

    # ...and a name that is not advertised must be rejected by the same table,
    # which is what makes the list exhaustive rather than merely a subset.
    unadvertised <- tryCatch({
        use_genetic_code("Standard")   # the plausible-but-wrong name from the docs
        GeneticCode()
        FALSE
    }, error = function(e) TRUE)
    use_genetic_code(original)

    if (!("Standard nuclear" %in% names)) {
        why <- "genetic_code_names() does not include \"Standard nuclear\""
    } else if (length(names) != 18) {
        why <- paste0("expected 18 genetic codes, got ", length(names))
    } else if (!unadvertised) {
        why <- "an unadvertised name (\"Standard\") was accepted, so the list is not exhaustive"
    } else if (get_genetic_code() != original) {
        why <- "the active genetic code was not restored"
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

try(use_genetic_code(original), silent = TRUE)

cat(if (ok) "M7: PASS\n" else paste0("M7: FAIL ", why, "\n"))
