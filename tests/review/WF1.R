# WF1 -- the reported equilibrium is the stationary distribution of the
# generator, under both fixation models (MutationSelection.cpp,
# CoEvolution.cpp, Palantir.cpp).
#
# transition() uses Kimura's (1962) eq. (10) fixation probability,
# u(s) = (1 - e^{-s}) / (1 - e^{-2Ns}), for which detailed balance gives
# u(s) / u(-s) = e^{(2N - 1)s}. equilibrium() and equilibrium_to_fitness()
# used e^{2Ns}, the exponent of the weak-selection form u(s) = s / (1 - e^{-2Ns})
# (eq. 11), so the reported equilibrium was not the stationary distribution of
# the generator: every frequency ratio was off by e^{f_j - f_i}, which is
# large at small N. The fix uses 2N - 1 under the default and adds
# fixation = "weak-selection", whose generator has the 2N exponent.
#
# Expected outcome:
#   PRE-FIX:  FAIL -- at N = 10 the default model violates detailed balance by
#             tens of percent, and the `fixation` argument does not exist.
#   POST-FIX: PASS --
#     (1) for N = 1, 2, 10, 100 and both fixation models (and the default),
#         the reported equilibrium satisfies detailed balance with the
#         transition matrix to 1e-10 relative, and pi' Q = 0;
#     (2) the same for CoEvolution at N = 10 with a general delta;
#     (3) under "weak-selection" every rate equals mu_ij S / (1 - e^{-S}) with
#         S = 2N (f_j - f_i), so models at N = 10 and N = 1000 with the same S
#         agree to 1e-12; under "diffusion" they differ;
#     (4) equilibrium_to_fitness() with the matching fixation model and
#         nucleotide equilibrium gives a model whose stationary amino-acid
#         marginal -- the null vector of the generator, solved for directly --
#         equals the target (TV below 1e-10);
#     (5) an unknown fixation name, and a population size of at most 1/2 under
#         the diffusion exponent 2N - 1, raise errors naming the argument.
suppressMessages(library(PalantiR))

result <- tryCatch({
    fail <- character(0)

    gc_codons <- GeneticCode()
    aa <- unname(gc_codons)
    aa_order <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
    aa_idx <- match(aa, aa_order)

    nuc_eq <- c(0.22, 0.28, 0.31, 0.19)  # T, C, A, G
    hky <- HasegawaKishinoYano(nuc_eq, transition_rate = 2.5, transversion_rate = 1)

    # Scaled amino-acid fitnesses 2N f spanning [-10, 10]
    S_aa <- seq(-10, 10, length.out = 20)[c(3, 17, 8, 12, 1, 20, 5, 14, 10, 19,
                                            2, 16, 7, 11, 4, 18, 6, 13, 9, 15)]

    # Largest relative violation of pi_i Q_ij = pi_j Q_ji over the moves
    db_violation <- function(m) {
        Q <- m$transition
        p <- as.numeric(m$equilibrium)
        flux <- p * Q                     # flux[i, j] = pi_i Q_ij
        off <- Q > 0 & row(Q) != col(Q)
        max(abs(flux[off] - t(flux)[off]) / flux[off])
    }
    # Largest |(pi' Q)_j| relative to the largest outflux pi_j |Q_jj|
    null_violation <- function(m) {
        p <- as.numeric(m$equilibrium)
        max(abs(crossprod(p, m$transition))) / max(abs(p * diag(m$transition)))
    }
    check_stationary <- function(m, label) {
        db <- db_violation(m)
        nv <- null_violation(m)
        if (!isTRUE(db <= 1e-10 && nv <= 1e-10)) {
            fail <<- c(fail, sprintf("%s: detailed-balance violation %.3g, pi'Q %.3g",
                                     label, db, nv))
        }
    }

    # --- (1) MutationSelection, default and both fixation models -------------
    for (N in c(1, 2, 10, 100)) {
        f <- S_aa / (2 * N)
        check_stationary(MutationSelection(N, 1, hky, f, "none"),
                         sprintf("MutationSelection default N=%d", N))
        for (fx in c("diffusion", "weak-selection")) {
            m <- MutationSelection(N, 1, hky, f, "none", fixation = fx)
            if (!identical(m$fixation, fx)) {
                fail <- c(fail, sprintf("model$fixation is %s, not %s",
                                        format(m$fixation), fx))
            }
            check_stationary(m, sprintf("MutationSelection %s N=%d", fx, N))
        }
    }

    # --- (2) CoEvolution with a general delta --------------------------------
    set.seed(11)
    delta <- matrix(round(runif(400, 0.5, 2), 3), 20, 20)
    N <- 10
    f1 <- S_aa / (4 * N)
    f2 <- rev(S_aa) / (4 * N)
    for (fx in c("diffusion", "weak-selection")) {
        co <- CoEvolution(N, 1, hky, f1, f2, delta, "none", fixation = fx)
        check_stationary(co, sprintf("CoEvolution %s N=%d", fx, N))
        rm(co); invisible(gc())
    }

    # --- (3) weak selection depends on S alone ----------------------------------
    neutral <- MutationSelection(10, 1, hky, rep(0, 20), "none")$transition  # mu_ij
    move <- neutral > 0 & row(neutral) != col(neutral)
    build <- function(N, fx) MutationSelection(N, 1, hky, S_aa / (2 * N), "none",
                                               fixation = fx)$transition
    S <- outer(S_aa[aa_idx], S_aa[aa_idx], function(a, b) b - a)
    h <- ifelse(S == 0, 1, S / -expm1(-S))
    weak10 <- build(10, "weak-selection")
    weak1000 <- build(1000, "weak-selection")
    rel_pred <- max(abs(weak10[move] - neutral[move] * h[move]) / (neutral[move] * h[move]))
    rel_n <- max(abs(weak10[move] - weak1000[move]) / weak1000[move])
    if (!isTRUE(rel_pred <= 1e-12)) {
        fail <- c(fail, sprintf("weak-selection rate != mu S/(1-e^-S) (rel %.3g)", rel_pred))
    }
    if (!isTRUE(rel_n <= 1e-12)) {
        fail <- c(fail, sprintf("weak-selection rates differ between N=10 and N=1000 at equal S (rel %.3g)",
                                rel_n))
    }
    diff10 <- build(10, "diffusion")
    diff1000 <- build(1000, "diffusion")
    if (!isTRUE(max(abs(diff10[move] - diff1000[move]) / diff1000[move]) > 1e-2)) {
        fail <- c(fail, "diffusion rates do not depend on N at equal S")
    }

    # --- (4) inversion round trip through the generator's null vector ----------
    target <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.03,
                0.04, 0.05, 0.06, 0.02, 0.03, 0.04, 0.05, 0.06, 0.04, 0.04)
    target <- target / sum(target)
    stationary_aa <- function(Q) {
        A <- t(Q)
        A[nrow(A), ] <- 1
        p <- solve(A, c(rep(0, nrow(A) - 1), 1))
        as.vector(tapply(p, factor(aa, levels = aa_order), sum))
    }
    for (N in c(10, 1000)) {
        for (fx in c("diffusion", "weak-selection")) {
            psi <- as.vector(equilibrium_to_fitness(target, N, nucleotide_equilibrium = nuc_eq,
                                                    fixation = fx))
            m <- MutationSelection(N, 1, hky, psi, "none", fixation = fx)
            tv <- 0.5 * sum(abs(stationary_aa(m$transition) - target))
            if (!isTRUE(tv <= 1e-10)) {
                fail <- c(fail, sprintf("round trip %s N=%d: TV %.3g from target", fx, N, tv))
            }
        }
    }

    # --- (5) argument checks ---------------------------------------------------
    errs <- function(expr, must_mention) {
        m <- tryCatch({ force(expr); NA_character_ },
                      error = function(e) conditionMessage(e))
        !is.na(m) && all(vapply(must_mention, function(w) grepl(w, m, fixed = TRUE), logical(1)))
    }
    if (!errs(MutationSelection(10, 1, hky, S_aa / 20, "none", fixation = "moran"),
              c("fixation", "moran"))) {
        fail <- c(fail, "an unknown fixation name is not rejected by MutationSelection")
    }
    if (!errs(equilibrium_to_fitness(target, 10, fixation = "moran"), c("fixation", "moran"))) {
        fail <- c(fail, "an unknown fixation name is not rejected by equilibrium_to_fitness")
    }
    if (!errs(equilibrium_to_fitness(target, 0.4), "population_size")) {
        fail <- c(fail, "population_size <= 1/2 is not rejected under the diffusion exponent 2N - 1")
    }
    if (errs(equilibrium_to_fitness(target, 0.4, fixation = "weak-selection"), "population_size")) {
        fail <- c(fail, "population_size = 0.4 is rejected under weak selection, whose exponent 2N is positive")
    }

    if (length(fail) == 0) "WF1: PASS" else paste0("WF1: FAIL ", paste(fail, collapse = "; "))
}, error = function(e) paste0("WF1: FAIL error: ", conditionMessage(e)))

cat(result, "\n", sep = "")
