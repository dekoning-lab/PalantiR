# E1 -- exact-method event positions through a stiff transient.
#
# When the entering forecast is far from the new mode's equilibrium, the scaled
# class flux can start orders of magnitude above its equilibrium value of ~1
# (634x measured on the disjoint-profile switch below). The exact rescaler's
# knot table used to march at a fixed intrinsic-time step, so the whole
# boundary layer fell inside one knot interval and the piecewise-linear
# inverse time change misplaced early events: 4.6x too many class events were
# reported by branch position 0.1, even though the events themselves, their
# total number, and the exit forecast were all correct.
#
# EXPECTED OUTCOME
#   Before the fix (adaptive knot refinement, 2026-08-21): FAIL. Cumulative
#     non-synonymous events per site by branch position u overshoot badly at
#     small u under rescale_method = "exact" (~0.46 at u = 0.1).
#   After the fix: PASS. Cumulative class events per site track u to within
#     0.04 at every checkpoint, under both rescale methods. The tolerance is
#     ~5 sampling SDs at the tightest checkpoint, and ~9x smaller than the
#     defect it guards against.

suppressMessages(library(PalantiR))

ID <- "E1"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))

result <- tryCatch({
    set_palantir_seed(20260821)
    data(psi_c50_1)
    hky <- HasegawaKishinoYano(rep(.25, 4))
    psiA <- equilibrium_to_fitness(as.numeric(psi_c50_1[1, ]), 1000)
    psiB <- equilibrium_to_fitness(as.numeric(psi_c50_1[2, ]), 1000)
    mA <- MutationSelection(1000, 1e-8, hky, psiA, scaling_type = "non-synonymous")
    mB <- MutationSelection(1000, 1e-8, hky, psiB, scaling_type = "non-synonymous")
    aa_of <- unname(GeneticCode())

    L <- 2; n_sites <- 1500
    treeT <- tempfile(fileext = ".newick"); writeLines(sprintf("(A:%f,B:0.01);", L), treeT)
    modeT <- tempfile(fileext = ".newick"); writeLines("(A:1,B:0);", modeT)
    pT <- Phylogeny(treeT); mT <- Phylogeny(modeT, type = "mode")

    for (method in c("segments", "exact")) {
        root <- sample_sequence(mA, n_sites)
        sim <- simulate_over_interval_phylogeny(pT, mT, list(mA, mB), root,
                                                start_mode = 0,
                                                rescale_method = method)
        nodeA <- sim$intervals$node[sim$intervals$state == 1]
        subA <- sim$substitutions[sim$substitutions$node %in% nodeA, ]
        ns_time <- subA$time[aa_of[subA$from + 1] != aa_of[subA$to + 1]]
        for (u in c(0.1, 0.25, 0.5, 1, 2)) {
            got <- sum(ns_time <= u) / n_sites
            if (abs(got - u) > 0.04) {
                fail(sprintf("(%s) cum non-syn/site %.3f at u=%.2f (expected %.2f +/- 0.04)",
                             method, got, u, u))
            }
        }
    }
    "PASS"
}, error = function(e) sprintf("FAIL (error: %s)", conditionMessage(e)))

done(sprintf("%s: %s", ID, result))
