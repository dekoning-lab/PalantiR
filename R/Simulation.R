plot.Simulation <- function(
    x,
    sim_type = x$type,
    plot_intervals = (x$type == "compound_codon"),
    ...) {
    simulation <- x
    # if simulation has multiple models
    PhyloPlot(simulation$phylogeny, simulation$substitutions, simulation$intervals,
              sim_type = sim_type,
              plot_intervals = plot_intervals, ...)
}

print.Simulation <- function(x, ...) {
    cat("Simulation of type", x$type, "with", ncol(x$alignment), "sites\n")
}
