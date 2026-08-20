# C7 -- seedable engine RNG.
#
# Before the fix: Palantir.hpp declared `static random_device rd; static
# mt19937 rng(rd());` at namespace scope in a header, so every translation unit
# had its own generator seeded from random_device. There was no seeding API at
# all, set.seed() had no effect, and no simulation was reproducible.
# EXPECTED against the OLD library: FAIL -- set_palantir_seed() does not exist
#   (`could not find function "set_palantir_seed"`).
# EXPECTED against the FIXED library: PASS -- two runs after the same
#   set_palantir_seed(1) are identical, and a different seed gives a different
#   result.

suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

tryCatch({
    tree <- Phylogeny(system.file("extdata", "primates.newick", package = "PalantiR"))
    hky  <- HasegawaKishinoYano(c(0.1, 0.2, 0.3, 0.4), transition_rate = 2)
    seq  <- Sequence(strrep("ACGT", 5), type = "nucleotide")

    run <- function(seed) {
        set_palantir_seed(seed)
        s <- simulate_over_phylogeny(tree, hky, seq, rate = 20)
        d <- s$substitutions
        # Waiting times come from Simulate.cpp, jump targets and the sampled
        # start states from Palantir.cpp: before the fix these were two
        # different engines, so both are part of the signature.
        list(n = nrow(d), time = d$time, to = d$to, from = d$from,
             aln = as.vector(unclass(s$alignment)),
             sampled = sample_sequence(hky, 25)$index)
    }

    a1 <- run(1)
    a2 <- run(1)
    b  <- run(2)

    if (!identical(a1, a2)) {
        why <- "two runs after set_palantir_seed(1) differ"
    } else if (a1$n == 0) {
        why <- "simulation produced no substitutions; signature is not informative"
    } else if (identical(a1, b)) {
        why <- "seed 1 and seed 2 gave identical results"
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "C7: PASS\n" else paste0("C7: FAIL ", why, "\n"))
