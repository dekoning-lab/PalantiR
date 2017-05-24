# Example of Co-Evolution simulation

# Load paramters
source("examples/common_parameters.R")

delta <- matrix(rep(1, 20 * 20), nrow = 20)

# Make Model
co <- CoEvolution(
    population_size = 1000,
    mutation_rate = 1e-8,
    nucleotide_mode = hky,
    fitness_1 = aa_psi[1,],
    fitness_2 = aa_psi[2,],
    delta = delta)

# sample sequence from model equilibrium
s <- sample_sequence(model = co, length = 100)

# simulate
sim <- simulate_over_phylogeny(
    phylogeny = mammals,
    model = co,
    sequence = s)
