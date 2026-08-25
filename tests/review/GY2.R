# GY2 -- exact scaling/omega identities and public-input validation for GY94.

suppressMessages(library(PalantiR))

ID <- "GY2"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

failures <- character(0)
check <- function(ok, message) if(!isTRUE(ok)) failures <<- c(failures, message)
errs <- function(expr, words) {
    message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
    !is.na(message) && all(vapply(words, grepl, logical(1), x = message, fixed = TRUE))
}

result <- try({
    original_code <- get_genetic_code()
    use_genetic_code("Standard nuclear")
    gc <- GeneticCode()
    codons <- names(gc)
    aa <- unname(gc)
    S <- length(codons)
    chars <- strsplit(codons, "")
    one_step <- outer(seq_len(S), seq_len(S), Vectorize(function(i, j) {
        i != j && sum(chars[[i]] != chars[[j]]) == 1L
    }))
    synonymous <- outer(aa, aa, "==") & one_step
    nonsynonymous <- !outer(aa, aa, "==") & one_step
    all_changes <- one_step
    class_rate <- function(model, mask) {
        sum(as.numeric(model$equilibrium) * rowSums(model$transition * mask))
    }

    pi <- F61(setNames(seq_len(S) + 10, codons))
    for(type in c("standard", "substitution", "synonymous", "non-synonymous")) {
        m <- GY94(.6, 2.7, pi, scaling_type = type)
        mask <- switch(type,
                       standard = all_changes,
                       substitution = all_changes,
                       synonymous = synonymous,
                       `non-synonymous` = nonsynonymous)
        check(abs(class_rate(m, mask) - 1) < 2e-12,
              sprintf("%s stationary class rate is not one", type))
    }
    m_standard <- GY94(.6, 2.7, pi, scaling_type = "standard")
    m_sub <- GY94(.6, 2.7, pi, scaling_type = "substitution")
    check(max(abs(m_standard$transition - m_sub$transition)) < 2e-12,
          "standard and substitution scaling aliases differ")

    lo <- .25
    hi <- 4
    syn_lo <- GY94(lo, 2.7, pi, scaling_type = "synonymous")
    syn_hi <- GY94(hi, 2.7, pi, scaling_type = "synonymous")
    check(abs(class_rate(syn_lo, synonymous) - 1) < 2e-12 &&
          abs(class_rate(syn_hi, synonymous) - 1) < 2e-12,
          "synonymous scaling does not pin synonymous flux")
    check(abs(class_rate(syn_hi, nonsynonymous) /
              class_rate(syn_lo, nonsynonymous) - hi / lo) < 2e-12,
          "nonsynonymous flux is not proportional to omega under synonymous scaling")

    ns_lo <- GY94(lo, 2.7, pi, scaling_type = "non-synonymous")
    ns_hi <- GY94(hi, 2.7, pi, scaling_type = "non-synonymous")
    check(abs(class_rate(ns_lo, nonsynonymous) - 1) < 2e-12 &&
          abs(class_rate(ns_hi, nonsynonymous) - 1) < 2e-12,
          "non-synonymous scaling does not pin nonsynonymous flux")
    check(abs(class_rate(ns_lo, synonymous) /
              class_rate(ns_hi, synonymous) - hi / lo) < 2e-12,
          "synonymous flux is not inversely proportional to omega under nonsynonymous scaling")
    ns_mask <- nonsynonymous
    check(max(abs((ns_lo$transition - ns_hi$transition)[ns_mask])) < 2e-12,
          "nonsynonymous Q entries changed with omega under nonsynonymous scaling")

    for(omega in c(lo, hi)) {
        m <- GY94(omega, 2.7, pi, scaling_type = "standard")
        check(abs(class_rate(m, all_changes) - 1) < 2e-12,
              "standard scaling does not pin total flux")
    }

    # Boundary and adapter validation: each failure must name its offending input.
    bad <- c(
        omega_zero = errs(GY94(omega = 0), "omega"),
        omega_inf = errs(GY94(omega = Inf), "omega"),
        kappa_zero = errs(GY94(kappa = 0), "kappa"),
        kappa_na = errs(GY94(kappa = NA_real_), "kappa"),
        scaling = errs(GY94(scaling_type = "totalish"), "scaling"),
        f1_length = errs(F1x4(c(.2, .3, .5)), "four"),
        f1_names = errs(F1x4(c(T=.25, C=.25, A=.25, U=.25)), c("T", "C", "A", "G")),
        f1_zero = errs(F1x4(c(.25, .25, .5, 0)), "positive"),
        f3_shape = errs(F3x4(matrix(1, 2, 6)), "3 x 4"),
        f3_nonfinite = errs(F3x4(matrix(c(rep(1, 11), Inf), 3, 4)), "positive"),
        f61_length = errs(F61(rep(1, S - 1L)), as.character(S)),
        f61_duplicate = errs(F61(setNames(rep(1, S), c(codons[-1], codons[2]))), "exactly once"),
        scheme_unknown = errs(codon_frequencies(c(.25,.25,.25,.25), "F2x4"), "Unknown"),
        scheme_disagrees = errs(codon_frequencies(FEqual(), "F61"), "disagrees"),
        site_both = errs(GY94SiteModel(list(a = m_standard), n_sites = 1,
                                       site_classes = "a"), "exactly one"),
        site_zero = errs(GY94SiteModel(list(a = m_standard), n_sites = 0), "positive"),
        site_unassigned = errs(GY94SiteModel(list(m_standard, m_sub),
                                             site_classes = c(1,1)), "Every model"),
        site_labels = errs(GY94SiteModel(list(m_standard, m_sub), n_sites = c(1,1),
                                        labels = c("x", "x")), "unique")
    )
    check(all(bad), paste("missing/unclear validation for", paste(names(bad)[!bad], collapse = ", ")))

    # A CodonFrequencies/model object cannot be reused after the global code changes.
    standard_pi <- FEqual()
    standard_model <- GY94(frequencies = standard_pi)
    use_genetic_code("Vertebrate mitochondrial")
    check(errs(codon_frequencies(standard_pi), c("Standard nuclear", "Vertebrate mitochondrial")),
          "CodonFrequencies genetic-code mismatch was accepted")
    check(errs(GY94BranchModel(list(standard_model),
                               mode_phylogeny = local({
                                   f <- tempfile(fileext = ".newick")
                                   writeLines("(A:0,B:0);", f)
                                   Phylogeny(f, type = "mode")
                               })), c("Standard nuclear", "Vertebrate mitochondrial")),
          "GY94 model genetic-code mismatch was accepted")
    use_genetic_code(original_code)

    if(length(failures)) fail(paste(failures, collapse = "; "))
    pass()
}, silent = TRUE)

try(use_genetic_code(if(exists("original_code")) original_code else "Standard nuclear"), silent = TRUE)
if(inherits(result, "try-error")) fail(gsub("[\r\n]+", " ", as.character(result)))
