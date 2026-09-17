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
    check(abs(flux(ms0, syn) - ms0$synonymous_opportunities) < 2e-12,
          "neutral mutation--selection synonymous rate does not equal N_S")

    # A dS-scaled mode change must run for ordinary branch time. Both rescaler
    # choices therefore take the same no-time-change path and consume the same
    # random stream.
    tf <- tempfile(fileext = ".newick")
    mf <- tempfile(fileext = ".newick")
    writeLines("(A:0.4,B:0.2);", tf)
    writeLines("(A:1,B:0);", mf)
    tree <- Phylogeny(tf)
    modes <- Phylogeny(mf, type = "mode")
    root <- sample_sequence(ms1, 40)
    set_palantir_seed(1701)
    seg <- simulate_over_interval_phylogeny(
        tree, modes, list(ms1, ms2), root, 0, rescale_method = "segments")
    set_palantir_seed(1701)
    ex <- simulate_over_interval_phylogeny(
        tree, modes, list(ms1, ms2), root, 0, rescale_method = "exact")
    check(identical(seg$alignment, ex$alignment) &&
          identical(seg$substitutions, ex$substitutions),
          "dS simulation was altered by the event-budget rescaler choice")

    if(length(failures)) fail(paste(failures, collapse = "; "))
    pass()
}, silent = TRUE)

if(inherits(result, "try-error")) fail(gsub("[\r\n]+", " ", as.character(result)))
