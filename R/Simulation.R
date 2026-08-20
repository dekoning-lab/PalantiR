# FIX (2026-08-20, M9): simulate_with_shared_time_heterogeneity() gives every
# site its own mode history, so its `intervals` member is a LIST of interval
# tables rather than one table (and since the M10 fix join() does the same). The
# widget's interval renderer filters a single flat table by node, so handed a
# list it matched nothing and drew a silently blank interval track. Resolve the
# list here rather than changing the js contract: a plot restricted to a single
# site gets that site's table, and anything wider gets no interval track plus a
# warning saying how to get one.
#
# `sites` are 0-based, as in the substitution table and the widget option; a
# joined simulation's interval list is named by the same index, so the lookup
# prefers the name and falls back to position for an unnamed per-site list.
.site_intervals <- function(intervals, site) {
    key <- as.character(site)
    if(!is.null(names(intervals)) && !is.na(match(key, names(intervals)))) {
        return(intervals[[key]])
    }
    index <- suppressWarnings(as.integer(site)) + 1L
    if(is.na(index) || index < 1 || index > length(intervals)) {
        warning("per-site intervals: site ", key, " has no interval table")
        return(NULL)
    }
    intervals[[index]]
}

.resolve_intervals <- function(intervals, options) {
    # one table (or none) already: nothing to resolve
    if(is.null(intervals) || is.data.frame(intervals) || !is.list(intervals)) {
        return(intervals)
    }

    sites <- if(is.na(match("sites", names(options)))) "all" else options[["sites"]]
    selected <- if(identical(sites, "all")) {
        seq_along(intervals) - 1L
    } else if(identical(sites, "none")) {
        integer(0)
    } else {
        sites
    }

    if(length(selected) == 0) {
        return(NULL)
    }
    if(length(selected) != 1) {
        warning("per-site intervals: pass sites=<one site> to display them")
        return(NULL)
    }
    .site_intervals(intervals, selected[[1]])
}

plot.Simulation <- function(
    x,
    sim_type = x$type,
    plot_intervals = (x$type == "compound_codon"),
    ...) {
    simulation <- x
    # if simulation has multiple models
    intervals <- .resolve_intervals(simulation$intervals, list(...))
    PhyloPlot(simulation$phylogeny, simulation$substitutions, intervals,
              sim_type = sim_type,
              plot_intervals = plot_intervals, ...)
}

print.Simulation <- function(x, ...) {
    cat("Simulation of type", x$type, "with", ncol(x$alignment), "sites\n")
}
