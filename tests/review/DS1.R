# DS1 -- dS scaling semantics for GY94 and mutation--selection models.
#
# The historical `synonymous` mode normalized one expected synonymous event
# per codon. That is not dS per synonymous opportunity and makes a supplied
# neutral nucleotide branch length too long by 1/N_S. The public
# `synonymous` name now selects the neutral dS gauge; the old event-count
# currency is retained as `synonymous-per-codon`.

suppressMessages(library(PalantiR))

ID <- "DS1"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

failures <- character(0)
check <- function(ok, message) if(!isTRUE(ok)) failures <<- c(failures, message)

result <- try({
    use_genetic_code("Standard nuclear")
    codons <- names(GeneticCode())
    aa <- unname(GeneticCode())
    chars <- strsplit(codons, "")
    S <- length(codons)
    one_step <- outer(seq_len(S), seq_len(S), Vectorize(function(i, j) {
        i != j && sum(chars[[i]] != chars[[j]]) == 1L
    }))
    syn <- outer(aa, aa, "==") & one_step
    flux <- function(model, mask) {
        sum(as.numeric(model$equilibrium) * rowSums(model$transition * mask))
    }

    nuc <- c(T = .22, C = .28, A = .31, G = .19)
    hky <- HasegawaKishinoYano(nuc, transition_rate = 2.5,
                               transversion_rate = 1)

    # GY94: omega-independent scale, neutral total rate three, and exact
    # conversion to the legacy per-codon synonymous-event gauge.
    pi <- F1x4(nuc)
    gy_a <- GY94(.2, 2.5, pi, scaling_type = "synonymous")
    gy_b <- GY94(4.0, 2.5, pi, scaling_type = "dS")
    gy_1 <- GY94(1.0, 2.5, pi, scaling_type = "ds")
    gy_c <- GY94(.2, 2.5, pi, scaling_type = "synonymous-per-codon")
    check(identical(gy_a$scaling_type, "dS") &&
          identical(gy_b$scaling_type, "dS") &&
          identical(gy_1$scaling_type, "dS"),
          "GY94 aliases did not canonicalize to dS")
    check(abs(gy_a$scaling - gy_b$scaling) < 1e-13,
          "GY94 dS denominator depends on omega")
    check(abs(flux(gy_1, one_step) - 3) < 2e-12,
          "neutral GY94 total rate is not three per codon")
    check(abs(sum(gy_1$equilibrium * gy_1$dS_outflux) - 1) < 2e-12,
          "GY94 neutral-reference dS clock does not have stationary rate one")
    check(abs(flux(gy_1, syn) - gy_1$synonymous_opportunities) < 2e-12,
          "GY94 synonymous rate does not equal N_S")
    check(abs(flux(gy_c, syn) - 1) < 2e-12,
          "legacy synonymous-per-codon rate is not one")
    check(max(abs(gy_a$transition -
                  gy_a$synonymous_opportunities * gy_c$transition)) < 2e-12,
          "GY94 dS and per-codon generators do not obey Q_dS = N_S Q_C")

    # Mutation--selection: the dS denominator comes from the mutation-only
    # reference and is shared by different amino-acid fitness profiles.
    N <- 200
    mu <- .01
    f1 <- seq(-.004, .004, length.out = 20)
    f2 <- rev(f1) + sin(seq_len(20)) * 1e-4
    ms1 <- MutationSelection(N, mu, hky, f1, "synonymous")
    ms2 <- MutationSelection(N, mu, hky, f2, "dS")
    ms0 <- MutationSelection(N, mu, hky, rep(0, 20), "dS")
    check(identical(ms1$scaling_type, "dS") &&
          identical(ms2$scaling_type, "dS"),
          "mutation--selection aliases did not canonicalize to dS")
    check(abs(ms1$scaling - ms2$scaling) < 1e-14,
          "mutation--selection dS denominator depends on the fitness profile")
    check(abs(flux(ms0, one_step) - 3) < 2e-12,
          "neutral mutation--selection total rate is not three per codon")
    check(abs(sum(ms0$equilibrium * ms0$dS_outflux) - 1) < 2e-12,
          "mutation--selection neutral-reference dS clock does not have rate one")
    check(abs(flux(ms0, syn) - ms0$synonymous_opportunities) < 2e-12,
          "neutral mutation--selection synonymous rate does not equal N_S")

    # A dS-scaled mode change uses the destination model's neutral-reference
    # per-nucleotide outflux as its clock functional. For neutral GY94 this is
    # exactly one third of the total event outflux, so a rescaled branch of
    # length one must deliver about three events per codon even when it enters
    # with a distribution far from the destination equilibrium.
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick")
        writeLines(txt, f)
        Phylogeny(f, type = type)
    }
    tree <- mk("(A:1);")
    modes <- mk("(A:1);", "mode")
    split <- floor(S / 2)
    w_a <- w_b <- rep(1, S)
    w_a[seq_len(split)] <- 1000
    w_b[seq_len(split)] <- .001
    gy_shift_a <- GY94(1, 4, F61(setNames(w_a, codons)),
                       scaling_type = "dS")
    gy_shift_b <- GY94(1, 4, F61(setNames(w_b, codons)),
                       scaling_type = "dS")
    check(sum(gy_shift_a$equilibrium * gy_shift_b$dS_outflux) < .3,
          "GY94 dS transient fixture does not require rescaling")
    root <- sample_sequence(gy_shift_a, 4000)
    for(method in c("segments", "exact")) {
        set_palantir_seed(1701)
        sim <- simulate_over_interval_phylogeny(
            tree, modes, list(gy_shift_a, gy_shift_b), root, 0,
            segment_length = .002, tolerance = 1e-9,
            rescale_method = method)
        events <- nrow(sim$substitutions) / 4000
        check(abs(events - 3) < .18,
              sprintf("%s dS rescaler delivered %.3f neutral events per codon, expected about 3",
                      method, events))
    }

    # Exercise the motivating mutation--selection cases directly: both the
    # amino-acid fitness profile and Ne change. The two numerical rescalers
    # should engage and approximate the same time change rather than silently
    # taking an ordinary branch-time path.
    aa_order <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
    p_a <- p_b <- setNames(rep(.05 / 19, 20), aa_order)
    p_a["W"] <- .95
    p_b["A"] <- .95
    n_a <- 500
    n_b <- 20000
    fit_a <- equilibrium_to_fitness(as.numeric(p_a), n_a,
                                    nucleotide_equilibrium = rep(.25, 4))
    fit_b <- equilibrium_to_fitness(as.numeric(p_b), n_b,
                                    nucleotide_equilibrium = rep(.25, 4))
    ms_a <- MutationSelection(n_a, mu, hky, fit_a, "dS")
    ms_b <- MutationSelection(n_b, mu, hky, fit_b, "dS")
    target_b <- sum(ms_b$equilibrium * ms_b$dS_outflux)
    entry_b <- sum(ms_a$equilibrium * ms_b$dS_outflux)
    check(abs(target_b - entry_b) > .1,
          "mutation--selection dS transient fixture does not require rescaling")
    root_ms <- sample_sequence(ms_a, 3000)
    set_palantir_seed(1702)
    seg <- simulate_over_interval_phylogeny(
        tree, modes, list(ms_a, ms_b), root_ms, 0,
        segment_length = .002, tolerance = 1e-9,
        rescale_method = "segments")
    set_palantir_seed(1702)
    ex <- simulate_over_interval_phylogeny(
        tree, modes, list(ms_a, ms_b), root_ms, 0,
        rescale_method = "exact")
    seg_events <- nrow(seg$substitutions) / 3000
    ex_events <- nrow(ex$substitutions) / 3000
    check(abs(seg_events - ex_events) < .15,
          sprintf("dS rescalers disagree after fitness/Ne shift: segments %.3f, exact %.3f",
                  seg_events, ex_events))
    check(!identical(seg$substitutions, ex$substitutions),
          "dS fitness/Ne shift bypassed both transient rescalers")

    if(length(failures)) fail(paste(failures, collapse = "; "))
    pass()
}, silent = TRUE)

if(inherits(result, "try-error")) fail(gsub("[\r\n]+", " ", as.character(result)))
