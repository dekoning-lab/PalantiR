# M4 -- Phylogeny() validates the newick file.
#
# Before the fix: Phylogeny() accepted anything ifstream could open. A
# directory, an empty file, arbitrary text and a file holding several trees all
# produced a "Phylogeny" object (the parser silently built a nonsense tree or
# used only the first tree).
# EXPECTED against the OLD library: FAIL -- the empty file, the garbage file,
#   the unbalanced tree and the two-tree file are all accepted.
# EXPECTED against the FIXED library: PASS -- each is rejected with a specific
#   message, and the trees shipped with the package still load.

suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

errs <- function(expr) {
    m <- tryCatch({ force(expr); NA_character_ },
                  error = function(e) conditionMessage(e))
    !is.na(m)
}

tryCatch({
    dir <- file.path(tempdir(), "M4")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    write_file <- function(name, text) {
        p <- file.path(dir, name); writeLines(text, p); p
    }

    empty      <- file.path(dir, "empty.newick"); file.create(empty)
    whitespace <- write_file("whitespace.newick", "   \n\t\n")
    garbage    <- write_file("garbage.newick", "this is not a tree at all")
    unbalanced <- write_file("unbalanced.newick", "((Human:1,Pan:1):1;")
    extra_open <- write_file("extra_close.newick", "(Human:1,Pan:1)):1;")
    two_trees  <- write_file("two.newick", "(Human:1,Pan:1);(Gorilla:1,Pongo:1);")
    trailing   <- write_file("trailing.newick", "(Human:1,Pan:1); junk")
    no_semi    <- write_file("no_semicolon.newick", "(Human:1,Pan:1)")

    checks <- c(
        "missing file"       = errs(Phylogeny(file.path(dir, "does_not_exist.newick"))),
        "directory"          = errs(Phylogeny(dir)),
        "empty file"         = errs(Phylogeny(empty)),
        "whitespace only"    = errs(Phylogeny(whitespace)),
        "garbage text"       = errs(Phylogeny(garbage)),
        "unclosed paren"     = errs(Phylogeny(unbalanced)),
        "unmatched close"    = errs(Phylogeny(extra_open)),
        "two trees"          = errs(Phylogeny(two_trees)),
        "trailing junk"      = errs(Phylogeny(trailing)),
        "no terminating ;"   = errs(Phylogeny(no_semi))
    )
    bad <- names(checks)[!checks]

    # Positive controls: every tree shipped with the package still loads, under
    # both types.
    shipped <- list.files(system.file("extdata", package = "PalantiR"), full.names = TRUE)
    shipped <- shipped[grepl("\\.(newick|tree)$", shipped)]
    loaded  <- vapply(shipped, function(f) !errs(Phylogeny(f)), logical(1))
    mode_ok <- !errs(Phylogeny(system.file("extdata", "primates_switch.newick",
                                           package = "PalantiR"), type = "mode"))

    if (length(bad) > 0) {
        why <- paste("no error for:", paste(bad, collapse = "; "))
    } else if (length(shipped) == 0) {
        why <- "found no shipped newick files to use as positive controls"
    } else if (!all(loaded)) {
        why <- paste("shipped tree rejected:",
                     paste(basename(shipped[!loaded]), collapse = ", "))
    } else if (!mode_ok) {
        why <- "Phylogeny(primates_switch.newick, type = \"mode\") was rejected"
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "M4: PASS\n" else paste0("M4: FAIL ", why, "\n"))
