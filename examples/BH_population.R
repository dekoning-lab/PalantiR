# Example of branch-heterogeneous simulation

# Load paramters
source("examples/common_parameters.R")

# Load branch labels
mammals_switches <- Phylogeny("data/mammals_switch.newick")

# Plot intervals
plot(mammals, mammals_switches)

# Create MS models
ms_1k <- MutationSelection(1000, 1e-8, hky, equilibrium_to_fitness(aa_pi[1, ], 1000))
ms_10k <- MutationSelection(10000, 1e-8, hky, equilibrium_to_fitness(aa_pi[1, ], 10000))

# Sample sequence
s <- sample_sequence(ms_1k, 100)

# Simulate
sim <- simulate_over_interval_phylogeny(mammals, mammals_switches, list(ms_1k, ms_10k), s, start_mode = 0)

