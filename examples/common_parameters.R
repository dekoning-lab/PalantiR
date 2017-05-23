#This script loads parameters common to all the examples

## Load mammalian phylogeny
mammals <- Phylogeny("data/mammals.newick")

## Create a "flat" nucleotide transition model
hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))

## Load amino acid equilibrium frequencies
aa_pi <- as.matrix(read.csv("data/amino_acid_equilibrium.csv"))

## Load amino acid fitness vectors
aa_psi <- as.matrix(read.csv("data/amino_acid_fitness_N_1000.csv"))
