# Example of Mutation Selection simulation

# Load paramters
source("examples/common_parameters.R")

# Make Model
ms <- MutationSelection(population_size = 1000, mutation_rate = 1e-8, nucleotide_model = hky, fitness = psi)

# sample sequence from model equilibrium
s <- sample_sequence(model = ms, length = 100)

# simulate
sim <- simulate_over_phylogeny(phylogeny = p, model = ms, sequence = s)

# examine alignment
plot(sim$alignment)


# examine substitutions
plot(sim, sites = 1:10)

# check specific substitution data
head(sim$substitutions)
