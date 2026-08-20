# MIN5 -- fixation_probability overflow guards (MutationSelection.cpp).
#
# For selection < ~-709.78, exp(-selection) overflows so numerator and
# denominator of (-expm1(-s))/(-expm1(-2Ns)) are both -Inf and the fixation
# probability is NaN, which propagates into the model's transition matrix.
# The fix returns 0 explicitly for selection <= -709 (where the unguarded
# formula already returns 0 whenever it is finite, because the denominator
# overflows first) and the 1 - e^{-s} asymptote for selection >= 709, leaving
# every normal-range value unchanged.
#
# Expected outcome:
#   PRE-FIX:  FAIL -- a model containing a fitness difference of -800 has NaN
#             transition entries.
#   POST-FIX: PASS --
#     (1) that model's transition matrix is all-finite, with rate exactly 0
#         into the lethal amino acid and rate 2N * mu (fixation probability
#         exactly 1) out of it;
#     (2) normal-range values are unchanged: implied fixation probabilities
#         over a grid of selection differences (|s| up to 25, |2Ns| up to 500)
#         match the unguarded formula to 1e-9 relative.
suppressMessages(library(PalantiR))

result <- tryCatch({
    fail <- character(0)

    gc_codons <- GeneticCode()
    S <- length(gc_codons)
    aa <- unname(gc_codons)
    aa_order <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

    nuc_eq <- c(0.22, 0.28, 0.31, 0.19)  # T, C, A, G
    hky <- HasegawaKishinoYano(nuc_eq, transition_rate = 2.5, transversion_rate = 1)

    # --- (1) extreme selection: -800 must give 0, not NaN --------------------
    N <- 100
    mu <- 0.01
    f_ext <- rep(0, 20)
    f_ext[match("L", aa_order)] <- -800    # below the -709.78 overflow point

    m0 <- MutationSelection(N, mu, hky, rep(0, 20), "none")  # neutral: Q0_ij = mu_ij
    mx <- MutationSelection(N, mu, hky, f_ext, "none")
    Q0 <- m0$transition
    Qx <- mx$transition
    ext <- which(aa == "L")
    nrm <- which(aa != "L")

    if (!isTRUE(all(is.finite(Qx)))) {
        fail <- c(fail, sprintf("%d non-finite transition entries at selection -800",
                                sum(!is.finite(Qx))))
    } else {
        into <- max(abs(Qx[nrm, ext]))
        if (!isTRUE(into == 0)) {
            fail <- c(fail, sprintf("rate into lethal amino acid %.6g != 0", into))
        }
        # out of the lethal class: fixation probability is exactly 1, so
        # Qx = 2N * mu_ij = 2N * Q0 there
        sel_out <- Q0[ext, nrm] > 0
        dev_out <- max(abs(Qx[ext, nrm][sel_out] / (2 * N * Q0[ext, nrm][sel_out]) - 1))
        if (!isTRUE(dev_out <= 1e-12)) {
            fail <- c(fail, sprintf("rate out of lethal amino acid != 2N * mu (rel dev %.6g)", dev_out))
        }
    }

    # --- (2) normal range unchanged ------------------------------------------
    N2 <- 10
    f_norm <- c(0, 0.05, -0.05, 0.2, -0.2, 0.5, -0.5, 1, -1, 2,
                -2, 3, -3, 5, -5, 8, -8, 12.5, -12.5, 0.1)
    m0b <- MutationSelection(N2, mu, hky, rep(0, 20), "none")
    mf <- MutationSelection(N2, mu, hky, f_norm, "none")

    # recover per-codon fitness (up to a constant) from the equilibria:
    # pi_f/pi_0 is proportional to exp(2N f), so differences of
    # log(pi_f/pi_0)/2N are the selection coefficients of each move
    d <- log(as.numeric(mf$equilibrium) / as.numeric(m0b$equilibrium)) / (2 * N2)
    sel <- outer(d, d, function(a, b) b - a)

    fix_ref <- function(N, s) ifelse(s == 0, 1 / (2 * N),
                                     (-expm1(-s)) / (-expm1(-2 * N * s)))
    pred <- 2 * N2 * m0b$transition * fix_ref(N2, sel)
    idx <- which(m0b$transition > 0 & row(pred) != col(pred))
    rel <- max(abs(mf$transition[idx] - pred[idx]) / pred[idx])
    if (!isTRUE(rel <= 1e-9)) {
        fail <- c(fail, sprintf("normal-range values changed (max rel dev %.6g)", rel))
    }

    if (length(fail) == 0) "MIN5: PASS" else paste0("MIN5: FAIL ", paste(fail, collapse = "; "))
}, error = function(e) paste0("MIN5: FAIL error: ", conditionMessage(e)))

cat(result, "\n", sep = "")
