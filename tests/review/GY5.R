# GY5 -- one common normalization across a GY94 site mixture, including one
# common factor per branch type in the Zhang--Nielsen--Yang model.

suppressMessages(library(PalantiR))

ID <- "GY5"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

failures <- character(0)
check <- function(ok, message) if(!isTRUE(ok)) failures <<- c(failures, message)
near <- function(x, y, tol = 2e-12) max(abs(as.numeric(x) - as.numeric(y))) < tol

result <- try({
    original_code <- get_genetic_code()
    use_genetic_code("Standard nuclear")
    gc <- GeneticCode()
    codons <- names(gc)
    amino_acids <- unname(gc)
    chars <- strsplit(codons, "")
    one_step <- outer(seq_along(codons), seq_along(codons),
                      Vectorize(function(i, j) {
                          i != j && sum(chars[[i]] != chars[[j]]) == 1L
                      }))
    nonsynonymous <- !outer(amino_acids, amino_acids, "==") & one_step
    all_changes <- one_step
    class_rate <- function(model, mask) {
        sum(as.numeric(model$equilibrium) * rowSums(model$transition * mask))
    }

    pi <- F1x4(c(T = .21, C = .29, A = .31, G = .19))
    counts <- c(slow = 3, fast = 1)
    weights <- counts / sum(counts)
    mixture <- GY94SiteModel(
        list(slow = GY94(.2, 2.3, pi, scaling_type = "standard"),
             fast = GY94(4.0, 2.3, pi, scaling_type = "standard")),
        n_sites = counts)
    rates <- vapply(mixture$models, class_rate, numeric(1), mask = all_changes)
    check(abs(sum(weights * rates) - 1) < 2e-12,
          "homogeneous mixture does not have weighted mean total rate one")
    check(abs(rates[[1]] - rates[[2]]) > .1,
          "homogeneous mixture classes were still normalized independently")
    check(length(unique(vapply(mixture$models, `[[`, numeric(1), "scaling"))) == 1L,
          "homogeneous classes do not share one scaling denominator")
    check(near(rates, mixture$stationary_class_rates[, 1]),
          "recorded class rates disagree with the normalized matrices")

    ns_mixture <- GY94SiteModel(
        list(low = GY94(.2, 2.3, pi, scaling_type = "non-synonymous"),
             high = GY94(4.0, 2.3, pi, scaling_type = "non-synonymous")),
        n_sites = c(1, 1))
    ns_rates <- vapply(ns_mixture$models, class_rate, numeric(1),
                       mask = nonsynonymous)
    check(abs(mean(ns_rates) - 1) < 2e-12,
          "non-synonymous mixture does not have weighted mean rate one")
    check(abs(ns_rates[[2]] / ns_rates[[1]] - 20) < 2e-12,
          "common non-synonymous scaling erased the omega rate difference")

    mode_path <- tempfile(fileext = ".newick")
    writeLines("(background:0,foreground:1);", mode_path)
    modes <- Phylogeny(mode_path, type = "mode")
    zny_counts <- c(2, 3, 4, 5)
    zny <- ZNYBranchSiteModel(
        modes, zny_counts, omega0 = .2, omega2 = 3.5, kappa = 2.3,
        frequencies = c(T = .15, C = .25, A = .35, G = .25),
        frequency_model = "F1x4", scaling_type = "standard")
    zny_weights <- zny_counts / sum(zny_counts)
    for(mode in 1:2) {
        models <- lapply(zny$models, function(specification) specification$models[[mode]])
        mode_rates <- vapply(models, class_rate, numeric(1), mask = all_changes)
        mode_scales <- vapply(models, `[[`, numeric(1), "scaling")
        check(abs(sum(zny_weights * mode_rates) - 1) < 2e-12,
              paste("ZNY weighted total rate is not one in mode", mode - 1L))
        check(length(unique(mode_scales)) == 1L,
              paste("ZNY classes do not share one scale in mode", mode - 1L))
        check(near(mode_scales[[1]], zny$mixture_scaling[[mode]]),
              paste("ZNY recorded mixture scale is wrong in mode", mode - 1L))
    }

    # Behavioural check of the exact transient rescaler. The foreground has a
    # different stationary distribution, so this exercises the class-specific
    # target rates used under one common mixture denominator.
    tree_path <- tempfile(fileext = ".newick")
    writeLines("(background:1,foreground:1);", tree_path)
    tree <- Phylogeny(tree_path)
    empirical <- ZNYBranchSiteModel(
        modes, rep(500, 4), omega0 = .2, omega2 = 3.5, kappa = 2.3,
        frequencies = c(T = .15, C = .25, A = .35, G = .25),
        frequency_model = "F1x4",
        foreground_frequencies = c(T = .40, C = .30, A = .20, G = .10),
        foreground_frequency_model = "F1x4", scaling_type = "standard")
    set_palantir_seed(20260825)
    sim <- simulate_gy94_site_model(tree, empirical, rescale_method = "exact")
    event_rate <- tabulate(sim$substitutions$node, nbins = 2L) / 2000
    check(all(abs(event_rate - 1) < .10),
          paste("simulated mixture event rates are not one:",
                paste(round(event_rate, 3), collapse = ", ")))
    check(!nrow(sim$substitutions) ||
          all(sim$substitutions$time >= 0 & sim$substitutions$time <= 1 + 1e-12),
          "mixture-scaled exact rescaling reported an event outside its branch")

    other_path <- tempfile(fileext = ".newick")
    writeLines("(background:1,foreground:0);", other_path)
    other_modes <- Phylogeny(other_path, type = "mode")
    incompatible <- tryCatch({
        GY94SiteModel(
            list(a = GY94BranchModel(list(GY94(.2), GY94(1)), modes),
                 b = GY94BranchModel(list(GY94(.2), GY94(1)), other_modes)),
            n_sites = c(1, 1))
        FALSE
    }, error = function(e) grepl("same mode phylogeny", conditionMessage(e), fixed = TRUE))
    check(incompatible, "incompatible site-class mode phylogenies were accepted")

    use_genetic_code(original_code)
    if(length(failures)) fail(paste(failures, collapse = "; "))
    pass()
}, silent = TRUE)

try(use_genetic_code(if(exists("original_code")) original_code else "Standard nuclear"), silent = TRUE)
if(inherits(result, "try-error")) fail(gsub("[\r\n]+", " ", as.character(result)))
