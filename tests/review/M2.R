# M2 -- degeneracy-aware equilibrium_to_fitness().
#
# Before the fix: equilibrium_to_fitness(pi, Ne) inverted pi -> psi ignoring the
# mutational weight of each amino acid's codons, so the amino-acid marginal of
# the resulting MutationSelection equilibrium was the codon-degeneracy-tilted
# transform of the target (a uniform target came out spread by exactly 6x).
# The fix adds an OPTIONAL `nucleotide_equilibrium` argument; the default
# (NULL) reproduces the old numbers bit for bit.
# EXPECTED against the OLD library: FAIL -- equilibrium_to_fitness() takes no
#   third argument ("unused argument").
# EXPECTED against the FIXED library: PASS -- (a) with nucleotide_equilibrium
#   supplied, the realised amino-acid marginal matches the target to a total
#   variation distance below 1e-6; (b) with it omitted, the returned fitness
#   vector is bit-identical to the values produced by the pre-fix library
#   (recorded below from rlib_patch); (c) the old tilt is still there when the
#   argument is omitted, so the two branches are genuinely different.

suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

# equilibrium_to_fitness(PI_REF, 1000) as produced by the pre-fix library
# (/Users/jasondk/bot-workspace/incidental-convergence/rlib_patch), printed with
# %.17g so that the comparison below is exact.
PI_REF <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.03,
            0.04, 0.05, 0.06, 0.02, 0.03, 0.04, 0.05, 0.06, 0.04, 0.04)
PI_REF <- PI_REF / sum(PI_REF)
PSI_OLD <- c(
    0.99919528104378297, 0.99939801359783698, 0.99954185463406298,
    0.9996534264097201,  0.99974458718811698, 0.99982166252803062,
    0.99988842822434287, 0.9999473197421711,  1,
    0.99939801359783698, 0.99954185463406298, 0.9996534264097201,
    0.99974458718811698, 0.99919528104378297, 0.99939801359783698,
    0.99954185463406298, 0.9996534264097201,  0.99974458718811698,
    0.99954185463406298, 0.99954185463406298)

tryCatch({
    Ne  <- 1000
    nuc <- rep(0.25, 4)
    hky <- HasegawaKishinoYano(nuc, transition_rate = 2)

    # Amino acids in PalantiR's index order, which is alphabetical: ACDEFGHIKLMNPQRSTVWY.
    aa_order <- sort(unique(as.vector(GeneticCode())))

    # Collapse a codon equilibrium to its amino-acid marginal, in aa_order.
    aa_marginal <- function(model) {
        aa <- as.vector(GeneticCode())
        m  <- tapply(as.vector(model$equilibrium), factor(aa, levels = aa_order), sum)
        as.vector(m)
    }
    tv <- function(a, b) 0.5 * sum(abs(a - b))

    # (a) round trip with the correction on: a uniform target and a ragged one.
    targets <- list(uniform = rep(1/20, 20), ragged = PI_REF)
    tvs <- vapply(targets, function(target) {
        psi <- as.vector(equilibrium_to_fitness(target, Ne, nucleotide_equilibrium = nuc))
        m   <- MutationSelection(Ne, 1, hky, psi, scaling_type = "non-synonymous")
        tv(aa_marginal(m), target)
    }, numeric(1))

    # (c) with the correction off, the uniform target is still tilted by the
    #     codon degeneracy (6-fold codon families vs 2-fold, etc.).
    psi_tilted <- as.vector(equilibrium_to_fitness(rep(1/20, 20), Ne))
    m_tilted   <- MutationSelection(Ne, 1, hky, psi_tilted, scaling_type = "non-synonymous")
    tilt       <- aa_marginal(m_tilted)
    tilt_ratio <- max(tilt) / min(tilt)

    # (b) the default path is unchanged, to the last bit.
    psi_default <- as.vector(equilibrium_to_fitness(PI_REF, Ne))

    if (any(tvs >= 1e-6)) {
        why <- paste0("round-trip TV distance too large: ",
                      paste(sprintf("%s=%.3g", names(tvs), tvs), collapse = ", "))
    } else if (!isTRUE(all.equal(psi_default, PSI_OLD, tolerance = 1e-15))) {
        why <- paste0("default (no nucleotide_equilibrium) path changed; max abs diff ",
                      sprintf("%.3g", max(abs(psi_default - PSI_OLD))))
    } else if (abs(tilt_ratio - 6) > 1e-6) {
        why <- paste0("expected the uncorrected uniform target to be spread 6x, got ",
                      sprintf("%.6f", tilt_ratio))
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "M2: PASS\n" else paste0("M2: FAIL ", why, "\n"))
