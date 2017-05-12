plot.Simulation <- function(simulation, ...) {
    # if simulation has multiple models
    PhyloPlot(simulation$phylogeny, simulation$substitutions, simulation$intervals,
              plot_intervals = length(simulation$models) > 1, ...)
}

print.Simulation <- function(simulation) {
    cat("Simulation of type", simulation$model$type, "with", ncol(simulation$alignment), "sites\n")
}
