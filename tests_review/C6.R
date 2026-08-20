# C6 -- MarkovModulated switching semantics (MarkovModulated.cpp).
#
# The old transition ignored the requested mode occupancy (the off-diagonal
# switching blocks used the bare symmetric exchangeability, making the
# realised occupancy uniform whatever was asked for) while
# switching_equilibrium instead DIVIDED the diagonal substitution blocks,
# inflating each mode's substitution rates by exactly 1/w. The fix builds the
# mode process as a proper GTR -- switching rate from mode i to mode j is
# S_ij * w_j -- and leaves each mode's substitution block as the component
# model's matrix unchanged.
#
# Expected outcome, with w = (0.8, 0.15, 0.05):
#   PRE-FIX:  FAIL -- solved occupancy is (1/3, 1/3, 1/3) and diagonal block 1
#             is inflated by 1/0.8 = 1.25.
#   POST-FIX: PASS --
#     (1) mode occupancy of the solved stationary distribution equals w to
#         1e-10;
#     (2) diagonal block 1's off-diagonal entries equal the component model's
#         transition entries exactly. (The design note says the block "equals
#         the component matrix exactly"; a proper generator must add the
#         mode-switch out-rate on the diagonal, so the diagonal is checked as
#         component diagonal minus sum_j S_1j w_j.)
#     (3) off-diagonal block (i, j) equals I * S_ij * w_j.
.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

result <- tryCatch({
    fail <- character(0)

    gc_codons <- GeneticCode()
    S <- length(gc_codons)

    nuc_eq <- c(0.22, 0.28, 0.31, 0.19)  # T, C, A, G
    hky <- HasegawaKishinoYano(nuc_eq, transition_rate = 2.5, transversion_rate = 1)

    N <- 100
    mu <- 0.01
    set.seed(1)
    f1 <- round(runif(20, -0.02, 0.02), 5)
    f2 <- round(runif(20, -0.02, 0.02), 5)
    f3 <- round(runif(20, -0.02, 0.02), 5)

    w <- c(0.8, 0.15, 0.05)
    S3 <- matrix(c(0.0, 1.3, 0.7,
                   1.3, 0.0, 2.1,
                   0.7, 2.1, 0.0), 3, 3)
    switching <- GeneralTimeReversible(w, S3)

    mods <- list(MutationSelection(N, mu, hky, f1, "synonymous"),
                 MutationSelection(N, mu, hky, f2, "synonymous"),
                 MutationSelection(N, mu, hky, f3, "synonymous"))

    mm <- MarkovModulatedMutationSelection(mods, switching)
    Q <- mm$transition
    pi <- as.numeric(mm$equilibrium)
    blk <- function(i, j) Q[((i - 1) * S + 1):(i * S), ((j - 1) * S + 1):(j * S)]

    # (1) solved stationary mode occupancy equals w
    occ <- sapply(1:3, function(k) sum(pi[((k - 1) * S + 1):(k * S)]))
    if (!isTRUE(max(abs(occ - w)) <= 1e-10)) {
        fail <- c(fail, sprintf("occupancy (%s) != w (max dev %.6g)",
                                paste(sprintf("%.4f", occ), collapse = ", "),
                                max(abs(occ - w))))
    }

    # (2) diagonal substitution blocks are the component matrices unchanged
    for (k in 1:3) {
        B <- blk(k, k)
        C <- mods[[k]]$transition
        off <- row(B) != col(B)
        dev_off <- max(abs(B[off] - C[off]))
        if (!isTRUE(dev_off <= 1e-14 * max(abs(C)))) {
            fail <- c(fail, sprintf("block (%d,%d) off-diagonal != component (max dev %.6g)",
                                    k, k, dev_off))
        }
        switch_out <- sum(S3[k, -k] * w[-k])
        dev_diag <- max(abs(diag(B) - (diag(C) - switch_out)))
        if (!isTRUE(dev_diag <= 1e-10 * max(abs(C)))) {
            fail <- c(fail, sprintf("block (%d,%d) diagonal != component diag - switch out-rate (max dev %.6g)",
                                    k, k, dev_diag))
        }
    }

    # (3) switching blocks: I * S_ij * w_j
    for (i in 1:3) for (j in 1:3) {
        if (i == j) next
        dev <- max(abs(blk(i, j) - diag(S) * (S3[i, j] * w[j])))
        if (!isTRUE(dev <= 1e-14)) {
            fail <- c(fail, sprintf("block (%d,%d) != I * S_ij * w_j (max dev %.6g)", i, j, dev))
        }
    }

    # sanity: still a generator
    if (!isTRUE(max(abs(rowSums(Q))) <= 1e-9)) {
        fail <- c(fail, "row sums of the compound generator are not zero")
    }

    if (length(fail) == 0) "C6: PASS" else paste0("C6: FAIL ", paste(fail, collapse = "; "))
}, error = function(e) paste0("C6: FAIL error: ", conditionMessage(e)))

cat(result, "\n", sep = "")
