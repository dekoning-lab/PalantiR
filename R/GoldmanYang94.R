# Goldman--Yang (GY94) codon models and discrete site classes.

.gy94_frequency_model <- function(x) {
    key <- tolower(gsub("[^[:alnum:]]", "", x))
    models <- c(fequal = "FEqual", f1x4 = "F1x4", f3x4 = "F3x4", f61 = "F61")
    value <- unname(models[key])
    if(length(value) != 1L || is.na(value)) {
        stop("Unknown codon-frequency model `", x,
             "`; choose FEqual, F1x4, F3x4, or F61")
    }
    value
}

.gy94_positive_frequencies <- function(x, argument) {
    if(!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
        stop("Argument `", argument, "` should contain strictly positive, finite frequencies")
    }
    x
}

.gy94_nucleotide_frequencies <- function(x, argument) {
    x <- .gy94_positive_frequencies(x, argument)
    if(length(x) != 4L) {
        stop("Argument `", argument, "` should have four nucleotide frequencies")
    }
    nucleotide_order <- c("T", "C", "A", "G")
    if(!is.null(names(x))) {
        if(anyDuplicated(names(x)) || !setequal(names(x), nucleotide_order)) {
            stop("Named nucleotide frequencies in `", argument,
                 "` should use exactly T, C, A, and G")
        }
        x <- x[nucleotide_order]
    } else {
        names(x) <- nucleotide_order
    }
    x / sum(x)
}

.new_codon_frequencies <- function(values, frequency_model, parameters = NULL) {
    codons <- names(GeneticCode())
    values <- as.numeric(values)
    names(values) <- codons
    structure(values,
              class = c("CodonFrequencies", "numeric"),
              frequency_model = frequency_model,
              frequency_parameters = parameters,
              genetic_code = get_genetic_code())
}

#' Equal frequencies over the sense codons in the active genetic code.
FEqual <- function() {
    codons <- names(GeneticCode())
    .new_codon_frequencies(rep(1 / length(codons), length(codons)), "FEqual")
}

#' Convert one set of nucleotide frequencies to conditional sense-codon frequencies.
F1x4 <- function(frequencies) {
    nucleotide <- .gy94_nucleotide_frequencies(frequencies, "frequencies")
    codons <- names(GeneticCode())
    values <- vapply(strsplit(codons, ""), function(codon) {
        prod(nucleotide[codon])
    }, numeric(1))
    values <- values / sum(values)
    .new_codon_frequencies(values, "F1x4", nucleotide)
}

#' Convert position-specific nucleotide frequencies to sense-codon frequencies.
F3x4 <- function(frequencies) {
    if(is.null(dim(frequencies)) && length(frequencies) == 12L) {
        frequencies <- matrix(frequencies, nrow = 3L, byrow = TRUE)
    }
    if(!is.matrix(frequencies) || !all(sort(dim(frequencies)) == c(3L, 4L))) {
        stop("Argument `frequencies` should be a 3 x 4 matrix (codon positions by T,C,A,G), ",
             "a 4 x 3 transpose, or a length-12 vector in position blocks")
    }

    nucleotide_order <- c("T", "C", "A", "G")
    if(nrow(frequencies) == 4L && ncol(frequencies) == 3L) {
        frequencies <- t(frequencies)
    }
    if(!is.null(colnames(frequencies))) {
        if(anyDuplicated(colnames(frequencies)) ||
           !setequal(colnames(frequencies), nucleotide_order)) {
            stop("The nucleotide dimension of `frequencies` should be named T, C, A, and G")
        }
        frequencies <- frequencies[, nucleotide_order, drop = FALSE]
    } else {
        colnames(frequencies) <- nucleotide_order
    }
    frequencies <- t(vapply(seq_len(3L), function(position) {
        .gy94_nucleotide_frequencies(frequencies[position, ], "frequencies")
    }, numeric(4L)))
    colnames(frequencies) <- nucleotide_order
    rownames(frequencies) <- paste0("position", seq_len(3L))

    codons <- names(GeneticCode())
    values <- vapply(strsplit(codons, ""), function(codon) {
        prod(frequencies[cbind(seq_len(3L), match(codon, nucleotide_order))])
    }, numeric(1))
    values <- values / sum(values)
    .new_codon_frequencies(values, "F3x4", frequencies)
}

#' Validate and normalize a full sense-codon frequency vector.
F61 <- function(frequencies) {
    frequencies <- .gy94_positive_frequencies(frequencies, "frequencies")
    codons <- names(GeneticCode())
    if(length(frequencies) != length(codons)) {
        stop("Argument `frequencies` should have ", length(codons),
             " entries under the active genetic code, not ", length(frequencies))
    }
    if(!is.null(names(frequencies))) {
        if(anyDuplicated(names(frequencies)) || !setequal(names(frequencies), codons)) {
            stop("Named F61 frequencies should contain every sense codon in the active genetic code exactly once")
        }
        frequencies <- frequencies[codons]
    }
    frequencies <- frequencies / sum(frequencies)
    .new_codon_frequencies(frequencies, "F61", frequencies)
}

#' Adapt any standard GY94 frequency specification to sense-codon frequencies.
codon_frequencies <- function(frequencies = NULL, frequency_model = NULL) {
    if(inherits(frequencies, "CodonFrequencies")) {
        if(!identical(attr(frequencies, "genetic_code"), get_genetic_code())) {
            stop("These codon frequencies were built under genetic code `",
                 attr(frequencies, "genetic_code"), "`, but the active code is `",
                 get_genetic_code(), "`")
        }
        if(!is.null(frequency_model) &&
           .gy94_frequency_model(frequency_model) != attr(frequencies, "frequency_model")) {
            stop("`frequency_model` disagrees with the supplied CodonFrequencies object")
        }
        return(frequencies)
    }

    if(is.null(frequency_model)) {
        if(is.null(frequencies)) {
            frequency_model <- "FEqual"
        } else if(is.matrix(frequencies) || length(frequencies) == 12L) {
            frequency_model <- "F3x4"
        } else if(length(frequencies) == 4L) {
            frequency_model <- "F1x4"
        } else if(length(frequencies) == length(GeneticCode())) {
            frequency_model <- "F61"
        } else {
            stop("Could not infer the codon-frequency model; specify `frequency_model`")
        }
    }
    frequency_model <- .gy94_frequency_model(frequency_model)
    switch(frequency_model,
           FEqual = {
               if(!is.null(frequencies) && length(frequencies) != 0L) {
                   stop("FEqual does not take frequency parameters")
               }
               FEqual()
           },
           F1x4 = F1x4(frequencies),
           F3x4 = F3x4(frequencies),
           F61 = F61(frequencies))
}

#' Construct a Goldman--Yang 1994 codon substitution model.
GoldmanYang94 <- function(omega = 1, kappa = 1, frequencies = NULL,
                          frequency_model = NULL, scaling_type = "standard") {
    pi <- codon_frequencies(frequencies, frequency_model)
    model <- .GoldmanYang94Cpp(
        equilibrium = unclass(pi),
        omega = omega,
        kappa = kappa,
        frequency_model = attr(pi, "frequency_model"),
        scaling_type = scaling_type)
    names(model$equilibrium) <- names(pi)
    names(model$codon_frequencies) <- names(pi)
    model$frequency_model <- attr(pi, "frequency_model")
    model$frequency_parameters <- attr(pi, "frequency_parameters")
    model
}

# Familiar short name used in the codon-model literature.
GY94 <- GoldmanYang94

.is_gy94 <- function(model) inherits(model, "GoldmanYang94") &&
    inherits(model, "SubstitutionModel")

.check_gy94_model <- function(model, argument = "model") {
    if(!.is_gy94(model)) {
        stop("Argument `", argument, "` should be a GoldmanYang94 substitution model")
    }
    if(!identical(model$genetic_code, get_genetic_code())) {
        stop("Argument `", argument, "` was built under genetic code `",
             model$genetic_code, "`, but the active code is `", get_genetic_code(), "`")
    }
    invisible(model)
}

.gy94_standalone_scaling <- function(model) {
    if(!is.null(model$standalone_scaling)) {
        return(as.numeric(model$standalone_scaling))
    }
    as.numeric(model$scaling)
}

.gy94_apply_mixture_scaling <- function(model, common_scaling) {
    standalone_scaling <- .gy94_standalone_scaling(model)
    current_scaling <- as.numeric(model$scaling)
    if(length(common_scaling) != 1L || !is.finite(common_scaling) ||
       common_scaling <= 0) {
        stop("The mixture-wide GY94 scaling factor should be a positive finite number")
    }

    # Recover the unnormalised GY94 generator, then apply the common
    # mixture-wide denominator. Multiplying a generator by a scalar does not
    # change its jump-destination probabilities, so `sampling` remains valid.
    model$transition <- model$transition * current_scaling / common_scaling
    model$standalone_scaling <- standalone_scaling
    model$scaling <- as.numeric(common_scaling)
    model$stationary_scaled_rate <- standalone_scaling / common_scaling
    model$scaling_scope <- "site_mixture"
    model
}

.gy94_common_site_scaling <- function(models, weights) {
    branch <- vapply(models, inherits, logical(1), "GY94BranchModel")
    reference <- if(any(branch)) models[[which(branch)[1L]]] else NULL
    n_modes <- if(is.null(reference)) 1L else length(reference$models)

    if(any(branch)) {
        reference_newick <- reference$mode_phylogeny$newick
        compatible <- vapply(models[branch], function(specification) {
            identical(specification$mode_phylogeny$newick, reference_newick) &&
                length(specification$models) == n_modes
        }, logical(1))
        if(!all(compatible)) {
            stop("All branch-heterogeneous classes in a GY94 site mixture should use ",
                 "the same mode phylogeny and number of branch models")
        }
    }

    underlying <- lapply(models, function(specification) {
        if(.is_gy94(specification)) list(specification) else specification$models
    })
    all_models <- unlist(underlying, recursive = FALSE)
    for(i in seq_along(all_models)) {
        .check_gy94_model(all_models[[i]], "site-class model")
    }
    scaling_types <- vapply(all_models, `[[`, character(1), "scaling_type")
    if(length(unique(scaling_types)) != 1L) {
        stop("All GY94 site classes should use the same `scaling_type`")
    }

    standalone <- do.call(rbind, lapply(models, function(specification) {
        values <- if(.is_gy94(specification)) {
            rep(.gy94_standalone_scaling(specification), n_modes)
        } else {
            vapply(specification$models, .gy94_standalone_scaling, numeric(1))
        }
        as.numeric(values)
    }))
    common <- colSums(standalone * weights)

    effective <- lapply(models, function(specification) {
        if(.is_gy94(specification)) {
            if(n_modes == 1L) {
                return(.gy94_apply_mixture_scaling(specification, common[[1L]]))
            }
            branch_models <- lapply(common, function(scale) {
                .gy94_apply_mixture_scaling(specification, scale)
            })
            GY94BranchModel(branch_models, reference$mode_phylogeny,
                            start_mode = reference$start_mode)
        } else {
            branch_models <- Map(.gy94_apply_mixture_scaling,
                                 specification$models, as.list(common))
            GY94BranchModel(branch_models, specification$mode_phylogeny,
                            start_mode = specification$start_mode)
        }
    })

    rownames(standalone) <- names(models)
    colnames(standalone) <- paste0("mode", seq_len(n_modes) - 1L)
    names(common) <- colnames(standalone)
    stationary_rates <- sweep(standalone, 2L, common, "/")
    list(models = effective,
         scaling_type = unique(scaling_types),
         common_scaling = common,
         stationary_rates = stationary_rates,
         mode_phylogeny = if(is.null(reference)) NULL else reference$mode_phylogeny)
}

#' Describe a branch-heterogeneous GY94 process for one site class.
GY94BranchModel <- function(models, mode_phylogeny, start_mode = 0L) {
    if(!is.list(models) || !length(models)) {
        stop("Argument `models` should contain at least one GY94 model")
    }
    for(i in seq_along(models)) .check_gy94_model(models[[i]], paste0("models[[", i, "]]"))
    scaling <- vapply(models, `[[`, character(1), "scaling_type")
    if(length(unique(scaling)) != 1L) {
        stop("All branch models should use the same `scaling_type`")
    }
    states <- vapply(models, `[[`, numeric(1), "n_states")
    if(length(unique(states)) != 1L) {
        stop("All branch models should use the same codon state space")
    }
    if(!inherits(mode_phylogeny, "Phylogeny") || !identical(mode_phylogeny$type, "mode")) {
        stop("Argument `mode_phylogeny` should be a mode Phylogeny")
    }
    if(length(start_mode) != 1L || is.na(start_mode) || start_mode != as.integer(start_mode) ||
       start_mode < 0L || start_mode >= length(models)) {
        stop("Argument `start_mode` should be a zero-based index into `models`")
    }
    structure(list(models = models, mode_phylogeny = mode_phylogeny,
                   start_mode = as.integer(start_mode)),
              class = "GY94BranchModel")
}

#' Define an arbitrary discrete collection of GY94 site classes.
GY94SiteModel <- function(models, n_sites = NULL, site_classes = NULL,
                          labels = names(models)) {
    if(!is.list(models) || !length(models)) {
        stop("Argument `models` should contain at least one site-class model")
    }
    valid <- vapply(models, function(model) {
        .is_gy94(model) || inherits(model, "GY94BranchModel")
    }, logical(1))
    if(!all(valid)) {
        stop("Each site class should contain a GY94 model or GY94BranchModel")
    }
    if(is.null(labels) || length(labels) == 0L) labels <- as.character(seq_along(models))
    if(length(labels) != length(models) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
        stop("Argument `labels` should contain one unique, non-empty label per model")
    }
    labels <- as.character(labels)

    if(xor(is.null(n_sites), is.null(site_classes)) == FALSE) {
        stop("Supply exactly one of `n_sites` (counts per class) or `site_classes` (one assignment per site)")
    }
    if(!is.null(n_sites)) {
        if(!is.numeric(n_sites) || length(n_sites) != length(models) || anyNA(n_sites) ||
           any(!is.finite(n_sites)) || any(n_sites != as.integer(n_sites)) || any(n_sites < 1L)) {
            stop("Argument `n_sites` should contain one positive whole-number count per site class")
        }
        assignment <- rep(seq_along(models), as.integer(n_sites))
    } else if(is.character(site_classes) || is.factor(site_classes)) {
        assignment <- match(as.character(site_classes), labels)
        if(anyNA(assignment)) stop("Every value in `site_classes` should match a class label")
    } else {
        if(!is.numeric(site_classes) || anyNA(site_classes) ||
           any(site_classes != as.integer(site_classes)) ||
           any(site_classes < 1L | site_classes > length(models))) {
            stop("Numeric `site_classes` should use one-based indices into `models`")
        }
        assignment <- as.integer(site_classes)
    }
    if(!length(assignment)) stop("At least one site should be assigned")
    if(!all(seq_along(models) %in% assignment)) {
        stop("Every model should be assigned to at least one site")
    }

    names(models) <- labels
    n_sites <- tabulate(assignment, nbins = length(models))
    weights <- n_sites / sum(n_sites)
    normalized <- .gy94_common_site_scaling(models, weights)

    structure(list(models = normalized$models,
                   labels = labels,
                   assignment = assignment,
                   n_sites = n_sites,
                   class_weights = setNames(as.numeric(weights), labels),
                   scaling_type = normalized$scaling_type,
                   mixture_scaling = normalized$common_scaling,
                   stationary_class_rates = normalized$stationary_rates,
                   mode_phylogeny = normalized$mode_phylogeny),
              class = "GY94SiteModel")
}

DiscreteGY94 <- GY94SiteModel

.gy94_root_model <- function(model) {
    if(.is_gy94(model)) return(model)
    model$models[[model$start_mode + 1L]]
}

.gy94_subsequence <- function(sequence, sites) {
    Sequence(paste(sequence$sequence[sites], collapse = ""), type = "codon")
}

.label_site_class_simulation <- function(sim, label, class_index) {
    n <- ncol(sim$alignment)
    metadata <- data.frame(
        site = seq_len(n) - 1L,
        site_class = rep(as.character(label), n),
        class_index = rep(as.integer(class_index), n),
        stringsAsFactors = FALSE)
    sim$site_classes <- metadata
    attr(sim$alignment, "site_classes") <- metadata
    sim$substitutions$site_class <- rep(as.character(label), nrow(sim$substitutions))
    sim$substitutions$class_index <- rep(as.integer(class_index), nrow(sim$substitutions))
    sim$site_class_models <- list()
    sim$site_class_models[[as.character(label)]] <- .sim_models(sim)
    sim
}

.restore_gy94_site_order <- function(sim, grouped_sites, site_model) {
    order_to_original <- match(seq_along(grouped_sites), grouped_sites)
    sim$alignment <- sim$alignment[, order_to_original, drop = FALSE]

    if(nrow(sim$substitutions)) {
        sim$substitutions$site <- grouped_sites[sim$substitutions$site + 1L] - 1L
        sim$substitutions <- sim$substitutions[order(sim$substitutions$site,
                                                     sim$substitutions$node,
                                                     sim$substitutions$time), , drop = FALSE]
        rownames(sim$substitutions) <- NULL
    }
    if(is.list(sim$intervals) && !is.data.frame(sim$intervals)) {
        sim$intervals <- sim$intervals[order_to_original]
        names(sim$intervals) <- as.character(seq_along(sim$intervals) - 1L)
    }

    metadata <- data.frame(
        site = seq_along(site_model$assignment) - 1L,
        site_class = site_model$labels[site_model$assignment],
        class_index = site_model$assignment - 1L,
        stringsAsFactors = FALSE)
    sim$site_classes <- metadata
    attr(sim$alignment, "site_classes") <- metadata
    sim$site_model <- site_model
    sim$site_class_models <- setNames(site_model$models, site_model$labels)
    sim
}

#' Simulate a GY94 alignment with arbitrary discrete site classes.
simulate_gy94_site_model <- function(
    phylogeny, site_model, sequence = NULL, rate = 1,
    segment_length = 0.001, tolerance = 0.001, rescale_method = "exact") {
    if(!inherits(site_model, "GY94SiteModel")) {
        stop("Argument `site_model` should be a GY94SiteModel")
    }
    for(i in seq_along(site_model$models)) {
        model <- site_model$models[[i]]
        if(.is_gy94(model)) {
            .check_gy94_model(model, paste0("site_model$models[[", i, "]]"))
        } else {
            for(j in seq_along(model$models)) {
                .check_gy94_model(model$models[[j]],
                                  paste0("site_model$models[[", i, "]]$models[[", j, "]]"))
            }
        }
    }
    total_sites <- length(site_model$assignment)
    if(!is.null(sequence)) {
        if(!inherits(sequence, "Sequence") || !identical(sequence$type, "codon") ||
           sequence$length != total_sites) {
            stop("Argument `sequence` should be a codon Sequence with one state per simulated site")
        }
    }

    simulations <- vector("list", length(site_model$models))
    grouped_sites <- integer(0)
    for(i in seq_along(site_model$models)) {
        sites <- which(site_model$assignment == i)
        grouped_sites <- c(grouped_sites, sites)
        specification <- site_model$models[[i]]
        root <- if(is.null(sequence)) {
            sample_sequence(.gy94_root_model(specification), length(sites))
        } else {
            .gy94_subsequence(sequence, sites)
        }
        simulations[[i]] <- if(.is_gy94(specification)) {
            simulate_over_phylogeny(phylogeny, specification, root, rate = rate)
        } else {
            scaling_targets <- vapply(specification$models, function(model) {
                if(is.null(model$stationary_scaled_rate)) 1 else model$stationary_scaled_rate
            }, numeric(1))
            simulate_over_interval_phylogeny(
                phylogeny = phylogeny,
                mode_phylogeny = specification$mode_phylogeny,
                models = specification$models,
                sequence = root,
                start_mode = specification$start_mode,
                rate = rate,
                segment_length = segment_length,
                tolerance = tolerance,
                rescale_method = rescale_method,
                scaling_targets = scaling_targets)
        }
        simulations[[i]] <- .label_site_class_simulation(
            simulations[[i]], site_model$labels[[i]], i - 1L)
    }

    combined <- do.call(join, simulations)
    .restore_gy94_site_order(combined, grouped_sites, site_model)
}

#' Build the four standard Zhang--Nielsen--Yang branch-site classes.
ZNYBranchSiteModel <- function(
    mode_phylogeny, site_counts, omega0 = 0.2, omega2 = 2, kappa = 2,
    frequencies = NULL, frequency_model = NULL,
    foreground_frequencies = frequencies,
    foreground_frequency_model = frequency_model,
    scaling_type = "standard", start_mode = 0L) {
    if(length(site_counts) != 4L) {
        stop("Argument `site_counts` should give counts for classes 0, 1, 2a, and 2b")
    }
    background_pi <- codon_frequencies(frequencies, frequency_model)
    foreground_pi <- codon_frequencies(foreground_frequencies,
                                        foreground_frequency_model)

    branch_class <- function(background_omega, foreground_omega) {
        GY94BranchModel(
            models = list(
                GoldmanYang94(background_omega, kappa, background_pi,
                              scaling_type = scaling_type),
                GoldmanYang94(foreground_omega, kappa, foreground_pi,
                              scaling_type = scaling_type)),
            mode_phylogeny = mode_phylogeny,
            start_mode = start_mode)
    }
    classes <- list(
        `0` = branch_class(omega0, omega0),
        `1` = branch_class(1, 1),
        `2a` = branch_class(omega0, omega2),
        `2b` = branch_class(1, omega2))
    GY94SiteModel(classes, n_sites = site_counts, labels = names(classes))
}

ZNYBranchSite <- ZNYBranchSiteModel
