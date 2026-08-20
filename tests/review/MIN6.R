# MIN6 -- sim$intervals is R's NULL, not the integer 0.
#
# Before the fix: simulate_over_phylogeny() built its result with
# `_["intervals"] = NULL`, where NULL is the C++ null pointer literal. Rcpp
# wraps that as the integer 0, so is.null(sim$intervals) was FALSE for every
# homogeneous simulation and code that branches on it took the wrong path.
# EXPECTED against the OLD library: FAIL -- sim$intervals is 0L.
# EXPECTED against the FIXED library: PASS -- sim$intervals is NULL and the
#   element is still present in the result list.

suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

tryCatch({
    tree <- Phylogeny(system.file("extdata", "primates.newick", package = "PalantiR"))
    hky  <- HasegawaKishinoYano(rep(0.25, 4), transition_rate = 2)
    seq  <- Sequence("ACGTACGT", type = "nucleotide")
    sim  <- simulate_over_phylogeny(tree, hky, seq, rate = 1)

    if (!("intervals" %in% names(sim))) {
        why <- "the result has no `intervals` element at all"
    } else if (!is.null(sim$intervals)) {
        why <- paste0("sim$intervals is not NULL: ",
                      paste(class(sim$intervals), collapse = "/"), " = ",
                      paste(utils::head(as.character(sim$intervals), 3), collapse = ","))
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "MIN6: PASS\n" else paste0("MIN6: FAIL ", why, "\n"))
