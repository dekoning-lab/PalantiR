# M3 -- population_size, mutation_rate and sample_sequence length are validated.
#
# Before the fix: these arguments were declared `unsigned long long` on the Rcpp
# boundary, so the R double was coerced before anything could inspect it.
# population_size = 0 gave an all-NaN model that printed normally, -1 wrapped to
# 1.8e19, 2.5 was silently truncated, and mutation_rate = 0 gave a degenerate
# rate matrix. sample_sequence(m, 0) returned an empty sequence.
# EXPECTED against the OLD library: FAIL -- none of these raise an error.
# EXPECTED against the FIXED library: PASS -- each raises an error naming the
#   argument, and valid values still work unchanged.

suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

errs <- function(expr, must_mention) {
    m <- tryCatch({ force(expr); NA_character_ },
                  error = function(e) conditionMessage(e))
    if (is.na(m)) return(FALSE)
    all(vapply(must_mention, function(w) grepl(w, m, fixed = TRUE), logical(1)))
}

tryCatch({
    hky <- HasegawaKishinoYano(rep(0.25, 4), transition_rate = 2)
    psi <- equilibrium_to_fitness(rep(1/20, 20), 1000)

    checks <- c(
        "Ne = 0"                 = errs(MutationSelection(0, 1, hky, psi), "population_size"),
        "Ne = -1"                = errs(MutationSelection(-1, 1, hky, psi), "population_size"),
        "Ne = 2.5"               = errs(MutationSelection(2.5, 1, hky, psi), "population_size"),
        "Ne = NA"                = errs(MutationSelection(NA_real_, 1, hky, psi), "population_size"),
        "mutation_rate = 0"      = errs(MutationSelection(1000, 0, hky, psi), "mutation_rate"),
        "mutation_rate = -1"     = errs(MutationSelection(1000, -1, hky, psi), "mutation_rate"),
        "CoEvolution Ne = 0"     = errs(CoEvolution(0, 1, hky, rep(1, 20), rep(1, 20),
                                                    matrix(1, 20, 20)), "population_size"),
        "e_to_f Ne = 0"          = errs(equilibrium_to_fitness(rep(1/20, 20), 0), "population_size")
    )

    m <- MutationSelection(1000, 1, hky, psi)
    checks <- c(checks,
        "sample_sequence length 0"  = errs(sample_sequence(m, 0), "length"),
        "sample_sequence length -1" = errs(sample_sequence(m, -1), "length"),
        "sample_sequence length 1.5" = errs(sample_sequence(m, 1.5), "length"))

    bad <- names(checks)[!checks]

    # Positive controls: valid values are unaffected.
    s <- sample_sequence(m, 7)
    if (length(bad) > 0) {
        why <- paste("no error for:", paste(bad, collapse = "; "))
    } else if (s$length != 7) {
        why <- paste0("positive control: sample_sequence(m, 7) has length ", s$length)
    } else if (!all(is.finite(m$equilibrium)) || abs(sum(m$equilibrium) - 1) > 1e-10) {
        why <- "positive control: MutationSelection(1000, 1, ...) equilibrium is not a distribution"
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "M3: PASS\n" else paste0("M3: FAIL ", why, "\n"))
