plot.Simulation <- function(
    simulation,
    sim_type = simulation$type,
    plot_intervals = (simulation$type == "compound_codon"),
    ...) {
    # if simulation has multiple models
    PhyloPlot(simulation$phylogeny, simulation$substitutions, simulation$intervals,
              sim_type = sim_type,
              plot_intervals = plot_intervals, ...)
}

print.Simulation <- function(simulation) {
    cat("Simulation of type", simulation$type, "with", ncol(simulation$alignment), "sites\n")
}
