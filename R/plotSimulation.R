plot.Simulation <- function(simulation) PhyloPlot(simulation$phylogeny, simulation$substitutions, simulation$intervals)

print.Simulation <- function(simulation) {
    cat("Simulation of type", simulation$model$type, "with", ncol(simulation$alignment), "sites\n")
}
