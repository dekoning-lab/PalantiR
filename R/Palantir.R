.globals <- new.env()

.globals$genetic_code_name <- "Standard nuclear"

# FIX (2026-08-20, M7): use_genetic_code() validated nothing. A misspelt name, a
# number, or NULL was stored happily and then poisoned every later call: the
# failure surfaced deep inside a model constructor as "Unknown genetic code" or,
# worse, as "Unexpected integer in codon sequence" from a simulation, with no
# hint that the setter was to blame. Validate here, against the names the C++
# GeneticCode_table actually holds (genetic_code_names(), so the two cannot
# drift), and name the alternatives in the message.
use_genetic_code <- function(name) {
    valid <- genetic_code_names()
    if(!is.character(name) || (length(name) != 1) || is.na(name)) {
        stop("Argument `name` should be a single genetic code name. Valid names are: ",
             paste(sort(valid), collapse = ", "))
    }
    if(is.na(match(name, valid))) {
        stop("Unknown genetic code \"", name, "\". Valid names are: ",
             paste(sort(valid), collapse = ", "))
    }
    .globals$genetic_code_name <- name
}

get_genetic_code <- function() .globals$genetic_code_name

.all_same <- function(x) all(x == x[[1]])
.get_member <- function(x, a) lapply(x, function(x) x[[a]])
.get_attr <- function(x, a) lapply(x, function(x) attr(x, a))

.pad_with_NA <- function(df, col_names) {
    for(name in col_names) {
        if(is.na(match(name, names(df)))) {
            df[[name]] <- NA
        }
    }
    df
}

# FIX (2026-08-20, M10): `site` in a substitution table is a 0-based index into
# that simulation's OWN alignment, so plain rbind() overlaid site 0 of every
# simulation on site 0 of the first. The alignments are cbind()ed in order, so
# each simulation's sites shift by the number of alignment columns contributed by
# the simulations before it.
.offset_sites <- function(df, offset) {
    if(!is.na(match("site", names(df)))) {
        df$site <- df$site + offset
    }
    df
}

# FIX (2026-08-20, M10/M9): join() used to concatenate the `intervals` members
# verbatim, giving one entry per SIMULATION while everything else in the joined
# result (alignment columns, substitution `site` indices) is per SITE. The joined
# structure is now one interval table per site of the joined alignment, named by
# the 0-based site index -- joined$intervals[["3"]] is the interval table of
# alignment column 4, the same index the `site` column and the widget's `sites`
# option use. The three sources normalize as follows:
#   - a homogeneous simulation has no intervals    -> NULL for each of its sites
#   - simulate_over_interval_phylogeny() and
#     simulate_with_shared_substitution_heterogeneity() share one table across
#     all sites                                    -> that table for each site
#   - simulate_with_shared_time_heterogeneity() is already per site
#                                                  -> element-wise
# If no simulation contributed a table the result is NULL, as for a homogeneous
# simulation.
.expand_intervals <- function(intervals, n_sites) {
    if(is.null(intervals) || !(is.data.frame(intervals) || is.list(intervals))) {
        # NULL, or the integer 0 that libraries built before the MIN6 fix used
        return(vector("list", n_sites))
    }
    if(is.data.frame(intervals)) {
        return(rep(list(intervals), n_sites))
    }
    if(length(intervals) != n_sites) {
        stop("A simulation has ", length(intervals), " interval tables for ",
             n_sites, " sites; expected one per site")
    }
    unname(as.list(intervals))
}

# A homogeneous simulation stores its single model under `model`; every
# heterogeneous one stores a list under `models`. .get_member(sims, "model")
# returned NULL for the latter, so joining interval simulations silently produced
# models = list(NULL, NULL). One entry per simulation either way.
.sim_models <- function(sim) if(is.null(sim$models)) sim$model else sim$models

.simulation_site_classes <- function(sim, n_sites, offset) {
    metadata <- sim$site_classes
    if(is.null(metadata)) {
        return(data.frame(
            site = seq_len(n_sites) - 1L + offset,
            site_class = rep(NA_character_, n_sites),
            class_index = rep(NA_integer_, n_sites),
            stringsAsFactors = FALSE))
    }
    if(!is.data.frame(metadata) || nrow(metadata) != n_sites ||
       !all(c("site", "site_class", "class_index") %in% names(metadata))) {
        stop("A simulation's `site_classes` should have one row per alignment column")
    }
    metadata$site <- metadata$site + offset
    metadata
}

# Join simulations
join <- function(...) {
    sims <- list(...)
    first <- sims[[1]]
    # check compatibiilty
    alignments <- .get_member(sims, "alignment")
    phylogenies <- .get_member(sims, "phylogeny")
    subs <- .get_member(sims, "substitutions")

    stopifnot(.all_same(.get_member(phylogenies, "newick")))
    stopifnot(.all_same(.get_attr(alignments, "type")))
    simulation_types <- unlist(.get_member(sims, "type"), use.names = FALSE)
    compatible_codon_types <- all(simulation_types %in% c("codon", "compound_codon"))
    if(!.all_same(simulation_types) && !compatible_codon_types) {
        stop("Simulations should have the same type (ordinary and interval codon simulations may be mixed)")
    }

    joined_alignment <- do.call(cbind, alignments)
    alignment_type <- attr(first$alignment, "type")
    attr(joined_alignment, "type") <- alignment_type
    class(joined_alignment) <- "Alignment"

    n_sites <- vapply(alignments, ncol, integer(1))
    offsets <- cumsum(c(0L, n_sites))[seq_along(sims)]

    col_names <- Reduce(union, .get_attr(subs, "names"))
    padded <- lapply(subs, .pad_with_NA, col_names)
    padded <- Map(.offset_sites, padded, offsets)
    joined_subs <- do.call(rbind, padded)

    joined_intervals <- unlist(
        Map(.expand_intervals, .get_member(sims, "intervals"), n_sites),
        recursive = FALSE)
    if(all(vapply(joined_intervals, is.null, logical(1)))) {
        joined_intervals <- NULL
    } else {
        names(joined_intervals) <- as.character(seq_along(joined_intervals) - 1L)
    }

    site_class_flags <- !vapply(.get_member(sims, "site_classes"), is.null, logical(1))
    has_site_classes <- any(site_class_flags)
    if(has_site_classes && !all(site_class_flags)) {
        stop("Either every joined simulation should carry `site_classes`, or none should")
    }
    joined_site_classes <- NULL
    if(has_site_classes) {
        joined_site_classes <- do.call(rbind, Map(
            .simulation_site_classes, sims, n_sites, offsets))
        rownames(joined_site_classes) <- NULL
        attr(joined_alignment, "site_classes") <- joined_site_classes
    }

    site_class_models <- unlist(.get_member(sims, "site_class_models"), recursive = FALSE)
    if(!length(site_class_models)) site_class_models <- NULL

    joined_type <- if(any(simulation_types == "compound_codon")) {
        "compound_codon"
    } else {
        first$type
    }

    joined_sim <- list(
        phylogeny = first$phylogeny,
        models = lapply(sims, .sim_models),
        substitutions = joined_subs,
        alignment = joined_alignment,
        intervals = joined_intervals,
        site_classes = joined_site_classes,
        site_class_models = site_class_models,
        type = joined_type
    )
    class(joined_sim) <- "Simulation"
    return(joined_sim)
}
