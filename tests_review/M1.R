# M1 -- CoEvolution simultaneous double substitutions (CoEvolution.cpp).
#
# The old transition gave a double substitution (both codons change at once)
# the rate 2N * mu_ik * mu_jl * fix -- the product of two mutation rates,
# dimensionally rate^2 -- so at mutation_rate = 1 ~54% of off-diagonal rate
# mass was double jumps and the NORMALISED model depended on the absolute
# mutation_rate. The fix sets the instantaneous double-substitution rate to
# zero (standard CTMC; compensatory pairs arise as fast sequential single
# moves through the delta coupling).
#
# Expected outcome:
#   PRE-FIX:  FAIL -- double-jump mass is nonzero, the normalised model
#             changes with mutation_rate, and the Kronecker identity is
#             violated.
#   POST-FIX: PASS --
#     (1) every transition entry where BOTH codons change is exactly 0;
#     (2) the normalised model at mutation_rate = 1 equals the one at
#         mutation_rate = 1e-6;
#     (3) with delta == 1 the generator equals the Kronecker sum of the two
#         independent MutationSelection generators divided by the mean of
#         their scaling factors (checks every entry, including the structural
#         zeros at double-jump positions).
.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

result <- tryCatch({
    fail <- character(0)
    tol <- 1e-8

    gc_codons <- GeneticCode()
    S <- length(gc_codons)

    nuc_eq <- c(0.22, 0.28, 0.31, 0.19)  # T, C, A, G
    hky <- HasegawaKishinoYano(nuc_eq, transition_rate = 2.5, transversion_rate = 1)

    N <- 100
    set.seed(1)
    f1 <- round(runif(20, -0.02, 0.02), 5)
    f2 <- round(runif(20, -0.02, 0.02), 5)
    delta_one <- matrix(1, 20, 20)

    # mutation_rate = 1: pre-fix, double jumps carry ~54% of the rate mass here
    co_a <- CoEvolution(N, 1, hky, f1, f2, delta_one, "synonymous")
    Qa <- co_a$transition
    pia <- as.numeric(co_a$equilibrium)
    rm(co_a); invisible(gc())

    # (1) zero instantaneous double-jump mass
    c1 <- rep(0:(S - 1), each = S)
    c2 <- rep(0:(S - 1), times = S)
    dbl <- outer(c1, c1, "!=") & outer(c2, c2, "!=")
    dbl_mass <- sum(abs(Qa[dbl]))
    if (!isTRUE(dbl_mass == 0)) {
        fail <- c(fail, sprintf("double-jump rate mass %.6g != 0", dbl_mass))
    }
    rm(dbl); invisible(gc())

    # (2) normalised model invariant to the absolute mutation_rate
    co_b <- CoEvolution(N, 1e-6, hky, f1, f2, delta_one, "synonymous")
    mu_dev <- max(abs(co_b$transition - Qa)) / max(abs(Qa))
    if (!isTRUE(mu_dev <= tol)) {
        fail <- c(fail, sprintf("normalised model depends on mutation_rate (rel dev %.6g)", mu_dev))
    }
    rm(co_b); invisible(gc())

    # (3) delta == 1: two independent MutationSelection sites (Kronecker sum)
    m1 <- MutationSelection(N, 1, hky, f1, "none")   # scaling "none" => unscaled Q
    m2 <- MutationSelection(N, 1, hky, f2, "none")
    r1 <- MutationSelection(N, 1, hky, f1, "synonymous")$scaling
    r2 <- MutationSelection(N, 1, hky, f2, "synonymous")$scaling
    I61 <- diag(S)
    Qpred <- (kronecker(m1$transition, I61) + kronecker(I61, m2$transition)) / ((r1 + r2) / 2)
    kron_dev <- max(abs(Qa - Qpred)) / max(abs(Qpred))
    if (!isTRUE(kron_dev <= tol)) {
        fail <- c(fail, sprintf("delta=1 generator != Kronecker sum of components (rel dev %.6g)", kron_dev))
    }
    pi_dev <- max(abs(pia - kronecker(as.numeric(m1$equilibrium), as.numeric(m2$equilibrium))))
    if (!isTRUE(pi_dev <= 1e-12)) {
        fail <- c(fail, sprintf("delta=1 equilibrium != product measure (max dev %.6g)", pi_dev))
    }

    if (length(fail) == 0) "M1: PASS" else paste0("M1: FAIL ", paste(fail, collapse = "; "))
}, error = function(e) paste0("M1: FAIL error: ", conditionMessage(e)))

cat(result, "\n", sep = "")
