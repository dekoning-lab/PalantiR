# M13 -- log-sum-exp offset and CoEvolution equilibrium renormalisation
# (MutationSelection.cpp, CoEvolution.cpp).
#
# The log-sum-exp offset was initialised to numeric_limits<double>::min(),
# the smallest POSITIVE double (2.2e-308), so an all-negative fitness vector
# never updated it: every term underflowed and the equilibrium came out
# all-NaN (MutationSelection) or non-finite (CoEvolution, which additionally
# never renormalised its equilibrium). The fix initialises the offset to
# -infinity and renormalises the CoEvolution equilibrium to sum to 1 exactly
# as MutationSelection does.
#
# Expected outcome:
#   PRE-FIX:  FAIL -- non-finite equilibria for all-negative fitness.
#   POST-FIX: PASS -- for all-negative fitness vectors both models yield
#             finite equilibria summing to 1, matching an R-side log-sum-exp
#             reference (and, for CoEvolution with delta == 1, the product
#             measure of the two site references).
.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

result <- tryCatch({
    fail <- character(0)

    gc_codons <- GeneticCode()
    S <- length(gc_codons)
    codons <- names(gc_codons)
    aa <- unname(gc_codons)
    aa_order <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
    aa_idx <- match(aa, aa_order)
    nuc_order <- c("T", "C", "A", "G")

    nuc_eq <- c(0.22, 0.28, 0.31, 0.19)  # T, C, A, G
    hky <- HasegawaKishinoYano(nuc_eq, transition_rate = 2.5, transversion_rate = 1)

    N <- 1000
    mu <- 1e-3
    # non-uniform, strictly negative fitness surfaces
    f1 <- seq(-600, -400, length.out = 20)
    f2 <- seq(-550, -420, length.out = 20)

    # R-side reference: pi_i proportional to exp(2N f_aa(i)) * prod_k nuc_eq
    cn <- t(sapply(codons, function(s) match(strsplit(s, "")[[1]], nuc_order)))
    lognuc <- rowSums(matrix(log(nuc_eq)[cn], S, 3))
    softmax <- function(lw) { w <- exp(lw - max(lw)); w / sum(w) }
    ref1 <- softmax(2 * N * f1[aa_idx] + lognuc)
    ref2 <- softmax(2 * N * f2[aa_idx] + lognuc)

    # MutationSelection
    ms <- MutationSelection(N, mu, hky, f1, "synonymous")
    eq <- as.numeric(ms$equilibrium)
    if (!isTRUE(all(is.finite(eq)))) {
        fail <- c(fail, "MutationSelection equilibrium not finite for all-negative fitness")
    } else {
        if (!isTRUE(abs(sum(eq) - 1) <= 1e-12)) {
            fail <- c(fail, sprintf("MutationSelection equilibrium sums to %.6g != 1", sum(eq)))
        }
        if (!isTRUE(max(abs(eq - ref1)) <= 1e-12)) {
            fail <- c(fail, sprintf("MutationSelection equilibrium off reference by %.6g",
                                    max(abs(eq - ref1))))
        }
    }

    # CoEvolution, delta == 1: equilibrium is the product measure ref1 x ref2
    delta_one <- matrix(1, 20, 20)
    co <- CoEvolution(N, mu, hky, f1, f2, delta_one, "synonymous")
    eqc <- as.numeric(co$equilibrium)
    if (!isTRUE(all(is.finite(eqc)))) {
        fail <- c(fail, "CoEvolution equilibrium not finite for all-negative fitness")
    } else {
        if (!isTRUE(abs(sum(eqc) - 1) <= 1e-12)) {
            fail <- c(fail, sprintf("CoEvolution equilibrium sums to %.6g != 1", sum(eqc)))
        }
        if (!isTRUE(max(abs(eqc - kronecker(ref1, ref2))) <= 1e-12)) {
            fail <- c(fail, sprintf("CoEvolution equilibrium off product-measure reference by %.6g",
                                    max(abs(eqc - kronecker(ref1, ref2)))))
        }
    }
    rm(co); invisible(gc())

    # CoEvolution with a non-trivial delta: still finite and normalised
    set.seed(2)
    delta_gen <- matrix(round(runif(400, 0.5, 2), 3), 20, 20)
    co2 <- CoEvolution(N, mu, hky, f1, f2, delta_gen, "synonymous")
    eqc2 <- as.numeric(co2$equilibrium)
    if (!isTRUE(all(is.finite(eqc2)) && abs(sum(eqc2) - 1) <= 1e-12)) {
        fail <- c(fail, "CoEvolution equilibrium with general delta not finite/normalised")
    }
    rm(co2); invisible(gc())

    if (length(fail) == 0) "M13: PASS" else paste0("M13: FAIL ", paste(fail, collapse = "; "))
}, error = function(e) paste0("M13: FAIL error: ", conditionMessage(e)))

cat(result, "\n", sep = "")
