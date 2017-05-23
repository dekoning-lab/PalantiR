# Example of population size switching simulation with Markov-Modulated Mutation Selection

# Load paramters
source("examples/common_parameters.R")

# Create Model

## parent "switching" model
gtr <- GeneralTimeReversible(
    equilibrium = flat(4),
    exchangeability = matrix(flat(16), nrow = 4))

# children MS models
pop_sizes <- c(100, 1000, 10000, 10000)
pop_models <- lapply(pop_sizes, function(pop_size) {
    pi <- aa_pi[5,]
    psi <- equilibrium_to_fitness(pi, pop_size)
    MutationSelection(pop_size, 1e-8, hky, psi)
})

# Markov-Modulated MS
mmms <- MarkovModulatedMutationSelection(pop_models, gtr)

# Sample sequence
s <- sample_sequence(mmms, 10)

# Simulate
sim <- simulate_over_phylogeny(p, mmms, s)
