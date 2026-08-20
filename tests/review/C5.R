# C5 -- CoEvolution::scaling rewrite (CoEvolution.cpp).
#
# The old scaling indexed the n_pairs x n_pairs pair transition matrix with
# single-codon indices (top-left 61x61 corner only), accumulated the first-site
# term inside the l loop and the second-site term inside the k loop (x61
# multiplicity), and classified the second site with _synonymous(i, k) instead
# of (j, l), so the time scale of every CoEvolution simulation was meaningless.
#
# Expected outcome:
#   PRE-FIX:  FAIL -- realised per-codon class event rate at equilibrium is far
#             from 1 for every scaling_type.
#   POST-FIX: PASS --
#     (1) for each scaling_type and for delta == 1 as well as a general delta,
#         the equilibrium-weighted class event rate per CODON site implied by
#         model$transition and model$equilibrium equals 1 (each single-site
#         move counts 1 if its changed codon is in the class; a pair state
#         covers two codon sites, hence the /2);
#     (2) with delta == 1 the model is two independent MutationSelection sites
#         and model$scaling equals the mean of the two component
#         MutationSelection scaling factors, whose own identity is known-exact.
suppressMessages(library(PalantiR))

result <- tryCatch({
    fail <- character(0)
    tol <- 1e-10

    gc_codons <- GeneticCode()
    S <- length(gc_codons)
    aa <- unname(gc_codons)

    nuc_eq <- c(0.22, 0.28, 0.31, 0.19)  # T, C, A, G
    hky <- HasegawaKishinoYano(nuc_eq, transition_rate = 2.5, transversion_rate = 1)

    N <- 100
    mu <- 0.01
    set.seed(1)
    # Weak selection (|2Nf| <= 4) so all three scaling types are
    # well-conditioned (strong selection drives the non-synonymous class
    # rate towards zero).
    f1 <- round(runif(20, -0.02, 0.02), 5)
    f2 <- round(runif(20, -0.02, 0.02), 5)
    delta_one <- matrix(1, 20, 20)
    set.seed(2)
    delta_gen <- matrix(round(runif(400, 0.5, 2), 3), 20, 20)

    # class matrix in codon space: cls[a, b] = TRUE if codon change a -> b is
    # an event of the requested class
    syn <- outer(aa, aa, "==")
    chg <- !diag(TRUE, S)
    class_of <- function(type) {
        if (type == "substitution") chg
        else if (type == "synonymous") chg & syn
        else chg & !syn
    }

    c1 <- rep(0:(S - 1), each = S)   # first-site codon of pair state p (0-based)
    c2 <- rep(0:(S - 1), times = S)  # second-site codon of pair state p

    # Expected class events per codon site per unit time at equilibrium:
    # sum_p pi_p sum_p' Q[p, p'] * (# changed codons of the move in the class),
    # divided by 2 codon sites per pair state. Aggregated per site through
    # rowsum() to avoid materialising 3721^2 helper matrices.
    realised_rate <- function(model, type) {
        A <- model$transition * as.numeric(model$equilibrium)  # A[p,p'] = pi_p Q[p,p']
        cls <- class_of(type)
        agg <- function(grp) t(rowsum(t(rowsum(A, grp)), grp)) # [a,b] = flow a -> b
        (sum(agg(c1) * cls) + sum(agg(c2) * cls)) / 2
    }

    for (type in c("substitution", "synonymous", "non-synonymous")) {
        ms1 <- MutationSelection(N, mu, hky, f1, type)
        ms2 <- MutationSelection(N, mu, hky, f2, type)

        co1 <- CoEvolution(N, mu, hky, f1, f2, delta_one, type)
        r <- realised_rate(co1, type)
        if (!isTRUE(abs(r - 1) <= tol)) {
            fail <- c(fail, sprintf("%s delta=1 realised per-codon rate %.6g != 1", type, r))
        }
        s_pred <- (ms1$scaling + ms2$scaling) / 2
        if (!isTRUE(abs(co1$scaling - s_pred) <= tol * abs(s_pred))) {
            fail <- c(fail, sprintf("%s delta=1 scaling %.6g != mean of component scalings %.6g",
                                    type, co1$scaling, s_pred))
        }
        rm(co1); invisible(gc())

        co2 <- CoEvolution(N, mu, hky, f1, f2, delta_gen, type)
        r2 <- realised_rate(co2, type)
        if (!isTRUE(abs(r2 - 1) <= tol)) {
            fail <- c(fail, sprintf("%s general-delta realised per-codon rate %.6g != 1", type, r2))
        }
        rm(co2); invisible(gc())
    }

    if (length(fail) == 0) "C5: PASS" else paste0("C5: FAIL ", paste(fail, collapse = "; "))
}, error = function(e) paste0("C5: FAIL error: ", conditionMessage(e)))

cat(result, "\n", sep = "")
