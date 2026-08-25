# GY1 -- deterministic GY94 codon-frequency adapters and rate matrix.
#
# This is deliberately an independent R construction of the Goldman--Yang
# generator.  It catches codon-order, target-frequency, transition and
# synonymous-classification errors without relying on simulation.

suppressMessages(library(PalantiR))

ID <- "GY1"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

failures <- character(0)
check <- function(ok, message) {
    if(!isTRUE(ok)) failures <<- c(failures, message)
}
near <- function(x, y, tol = 2e-12) {
    isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol,
                     check.attributes = FALSE))
}

result <- try({
    original_code <- get_genetic_code()
    use_genetic_code("Standard nuclear")
    gc <- GeneticCode()
    codons <- names(gc)
    aa <- unname(gc)
    S <- length(codons)

    # The four adapters coincide in their common uniform special case.
    equal <- FEqual()
    one <- F1x4(c(T = .25, C = .25, A = .25, G = .25))
    three <- F3x4(matrix(.25, 3, 4,
                         dimnames = list(NULL, c("T", "C", "A", "G"))))
    full <- F61(setNames(rep(1, S), rev(codons)))
    check(near(equal, rep(1 / S, S)), "FEqual is not uniform over sense codons")
    check(near(equal, one) && near(equal, three) && near(equal, full),
          "uniform FEqual/F1x4/F3x4/F61 adapters disagree")
    check(identical(names(equal), codons), "frequency vectors are not in GeneticCode() order")

    # F1x4 is the product distribution conditioned on being a sense codon.
    f1 <- c(T = .10, C = .20, A = .30, G = .40)
    expected1 <- vapply(strsplit(codons, ""), function(x) prod(f1[x]), numeric(1))
    expected1 <- expected1 / sum(expected1)
    pi1 <- F1x4(f1[c("G", "A", "T", "C")]) # named input may be permuted
    check(near(pi1, expected1), "F1x4 does not equal the conditioned codon-product distribution")
    check(near(codon_frequencies(f1), pi1), "codon_frequencies() did not infer F1x4")

    # Asymmetric rows expose a position transpose or a T,C,A,G ordering error.
    f3 <- rbind(
        position1 = c(T = .10, C = .20, A = .30, G = .40),
        position2 = c(T = .40, C = .30, A = .20, G = .10),
        position3 = c(T = .15, C = .25, A = .35, G = .25))
    expected3 <- vapply(strsplit(codons, ""), function(x) {
        prod(f3[cbind(1:3, match(x, colnames(f3)))])
    }, numeric(1))
    expected3 <- expected3 / sum(expected3)
    pi3 <- F3x4(f3)
    check(near(pi3, expected3), "F3x4 position-specific product is wrong")
    check(near(F3x4(t(f3)), pi3), "4 x 3 F3x4 transpose is not handled correctly")
    check(near(F3x4(as.vector(t(f3))), pi3), "length-12 F3x4 position blocks are wrong")
    check(near(codon_frequencies(f3), pi3), "codon_frequencies() did not infer F3x4")

    # Named F61 must be reordered, not merely have its names replaced.
    weights <- setNames(seq_len(S), codons)
    pi61 <- F61(weights[sample(seq_len(S))])
    check(near(pi61, weights / sum(weights)), "named F61 was not reordered to GeneticCode() order")
    check(identical(attr(pi61, "frequency_model"), "F61"), "F61 metadata is missing")

    # Independent reference Q for a deliberately asymmetric F61 distribution.
    omega <- 0.37
    kappa <- 3.4
    model <- GY94(omega = omega, kappa = kappa, frequencies = pi61,
                  scaling_type = "standard")
    chars <- strsplit(codons, "")
    is_transition <- function(a, b) paste0(a, b) %in% c("TC", "CT", "AG", "GA")
    raw <- matrix(0, S, S)
    for(i in seq_len(S)) for(j in seq_len(S)) if(i != j) {
        changed <- which(chars[[i]] != chars[[j]])
        if(length(changed) == 1L) {
            multiplier <- 1
            if(is_transition(chars[[i]][changed], chars[[j]][changed])) multiplier <- multiplier * kappa
            if(aa[i] != aa[j]) multiplier <- multiplier * omega
            raw[i, j] <- as.numeric(pi61[j]) * multiplier
        }
    }
    diag(raw) <- -rowSums(raw)
    rho <- sum(as.numeric(pi61) * -diag(raw))
    reference <- raw / rho

    check(max(abs(model$transition - reference)) < 2e-12,
          sprintf("GY94 matrix differs from independent reference (max %.3g)",
                  max(abs(model$transition - reference))))
    check(max(abs(rowSums(model$transition))) < 2e-12, "GY94 rows do not sum to zero")
    check(max(abs(as.numeric(pi61 %*% model$transition))) < 2e-12,
          "supplied codon frequencies are not stationary")
    detailed <- sweep(model$transition, 1, as.numeric(pi61), "*")
    check(max(abs(detailed - t(detailed))) < 2e-12, "GY94 detailed balance is violated")
    check(identical(model$scaling_type, "substitution"),
          "user-facing standard scaling was not canonicalized to substitution")
    check(identical(model$model, "GY94") && inherits(model, "GoldmanYang94") &&
          inherits(model, "SubstitutionModel"), "GY94 model identity/class metadata is wrong")
    check(near(GoldmanYang94(omega, kappa, pi61)$transition, model$transition),
          "GY94 and GoldmanYang94 aliases disagree")

    use_genetic_code(original_code)
    if(length(failures)) fail(paste(failures, collapse = "; "))
    pass()
}, silent = TRUE)

try(use_genetic_code(if(exists("original_code")) original_code else "Standard nuclear"), silent = TRUE)
if(inherits(result, "try-error")) fail(gsub("[\r\n]+", " ", as.character(result)))
