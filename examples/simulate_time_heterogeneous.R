#!/usr/bin/env Rscript
#
# simulate_time_heterogeneous.R
#
# Command-line driver for a time-heterogeneous mutation-selection simulation:
# the phylogeny is paired with a mode tree that assigns each branch a model
# index, and each mode runs a MutationSelection model with its own effective
# population size. The default inputs reproduce the documentation example: a
# mammalian phylogeny in which the primate clade evolves under a larger
# population size than the rest of the tree.
#
# Usage:
#   Rscript simulate_time_heterogeneous.R [--option=value ...]
#
# Options (defaults in parentheses):
#   --sites=N              number of alignment sites to simulate (100)
#   --population-sizes=A,B comma-separated effective population sizes, one per
#                          mode; branches labelled mode m in the mode tree use
#                          the (m+1)-th value (5000,30000)
#   --mutation-rate=X      per-generation mutation rate (1e-8)
#   --tree=FILE            newick tree with branch lengths
#                          (mammals.newick shipped with the package)
#   --mode-tree=FILE       newick tree with mode indices in place of branch
#                          lengths (mammals_switch.newick, which labels the
#                          primate clade mode 1 and all other branches mode 0)
#   --start-mode=M         mode at the root (0)
#   --profile=P            row of the packaged psi_c50_1 profile set to use as
#                          the amino-acid equilibrium for all sites (1)
#   --seed=N               random seed (20260820)
#   --out=DIR              output directory (time_heterogeneous_output)
#
# The amino-acid profile is converted to fitness at the root-mode population
# size, so the root sequence is sampled from the stationary distribution of
# the model the simulation starts in.
#
# Outputs, written to --out:
#   alignment_codon.fasta       simulated alignment, codon sequences
#   alignment_amino_acid.fasta  the same alignment, translated
#   substitutions.tsv           full substitution history (one row per event)
#   branch_modes.tsv            mode assigned to each branch (node index)
#   substitution_history.html   interactive substitution-history plot
#   run_info.txt                parameters, seed, versions, output summary
#
# Example:
#   Rscript simulate_time_heterogeneous.R --population-sizes=5000,8000 --sites=300

suppressMessages(library(PalantiR))

# ---------------------------------------------------------------------------
# Command-line options
# ---------------------------------------------------------------------------

defaults <- list(
    sites            = "100",
    population_sizes = "5000,30000",
    mutation_rate    = "1e-8",
    tree             = system.file("extdata", "mammals.newick", package = "PalantiR"),
    mode_tree        = system.file("extdata", "mammals_switch.newick", package = "PalantiR"),
    start_mode       = "0",
    profile          = "1",
    seed             = "20260820",
    out              = "time_heterogeneous_output")

argv <- commandArgs(trailingOnly = TRUE)
if (any(argv %in% c("-h", "--help"))) {
    cat("Options (--option=value):\n")
    for (k in names(defaults))
        cat(sprintf("  --%-18s (default: %s)\n", gsub("_", "-", k), defaults[[k]]))
    quit(status = 0)
}
opts <- defaults
for (a in argv) {
    if (!grepl("^--[a-z-]+=.+", a)) stop("unrecognized argument: ", a, " (see --help)")
    key <- gsub("-", "_", sub("=.*", "", sub("^--", "", a)))
    if (!key %in% names(opts)) stop("unknown option in: ", a, " (see --help)")
    opts[[key]] <- sub("^[^=]*=", "", a)
}

n_sites    <- as.integer(opts$sites)
pop_sizes  <- as.numeric(strsplit(opts$population_sizes, ",")[[1]])
mu         <- as.numeric(opts$mutation_rate)
start_mode <- as.integer(opts$start_mode)
profile    <- as.integer(opts$profile)
seed       <- as.integer(opts$seed)
stopifnot(n_sites > 0, length(pop_sizes) >= 2, all(pop_sizes > 0), mu > 0,
          start_mode >= 0, start_mode < length(pop_sizes))

set_palantir_seed(seed)

# ---------------------------------------------------------------------------
# Inputs: trees and fitness profile
# ---------------------------------------------------------------------------

phylogeny <- Phylogeny(opts$tree)
mode_tree <- Phylogeny(opts$mode_tree, type = "mode")

data(psi_c50_1, package = "PalantiR", envir = environment())
if (profile < 1 || profile > nrow(psi_c50_1))
    stop("--profile must be between 1 and ", nrow(psi_c50_1))

# One fitness vector for all sites and all modes, converted at the root-mode
# population size; the modes differ only in effective population size.
psi <- equilibrium_to_fitness(as.numeric(psi_c50_1[profile, ]),
                              pop_sizes[start_mode + 1])

# ---------------------------------------------------------------------------
# Models: one MutationSelection model per mode
# ---------------------------------------------------------------------------

nucleotide_model <- HasegawaKishinoYano(equilibrium = rep(0.25, 4))

models <- lapply(pop_sizes, function(N)
    MutationSelection(
        population_size  = N,
        mutation_rate    = mu,
        nucleotide_model = nucleotide_model,
        fitness          = psi))

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

message(sprintf("simulating %d sites, modes at N = %s, root mode %d ...",
                n_sites, paste(pop_sizes, collapse = ", "), start_mode))

root <- sample_sequence(models[[start_mode + 1]], n_sites)

sim <- simulate_over_interval_phylogeny(
    phylogeny      = phylogeny,
    mode_phylogeny = mode_tree,
    models         = models,
    sequence       = root,
    start_mode     = start_mode)

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

dir.create(opts$out, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(opts$out, f)

as.fasta(sim$alignment, file = out("alignment_codon.fasta"))

# as_amino_acid() translates elementwise but returns a plain vector; assigning
# through [] keeps the alignment's dimensions and taxon names.
aa_alignment <- sim$alignment
aa_alignment[] <- as_amino_acid(sim$alignment)
as.fasta(aa_alignment, file = out("alignment_amino_acid.fasta"))

write.table(sim$substitutions, out("substitutions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# sim$intervals holds one row per branch: node index, mode ("state"), and the
# span of the branch the mode covers.
write.table(sim$intervals, out("branch_modes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

widget_ok <- tryCatch({
    path <- normalizePath(out("substitution_history.html"), mustWork = FALSE)
    htmlwidgets::saveWidget(plot(sim), path)
    # saveWidget embeds all dependencies but can leave its staging directory
    unlink(sub("\\.html$", "_files", path), recursive = TRUE)
    TRUE
}, error = function(e) {
    message("substitution_history.html not written (", conditionMessage(e), ")")
    FALSE
})

writeLines(c(
    "PalantiR time-heterogeneous simulation",
    sprintf("date:              %s", format(Sys.time())),
    sprintf("PalantiR version:  %s", as.character(packageVersion("PalantiR"))),
    sprintf("R version:         %s", R.version.string),
    sprintf("tree:              %s", opts$tree),
    sprintf("mode tree:         %s", opts$mode_tree),
    sprintf("population sizes:  %s (mode 0..%d)", paste(pop_sizes, collapse = ", "),
            length(pop_sizes) - 1),
    sprintf("start mode:        %d", start_mode),
    sprintf("fitness profile:   psi_c50_1 row %d, converted at N = %g",
            profile, pop_sizes[start_mode + 1]),
    sprintf("sites:             %d", n_sites),
    sprintf("mutation rate:     %g", mu),
    sprintf("seed:              %d", seed),
    sprintf("substitutions:     %d events (%d non-synonymous)",
            nrow(sim$substitutions), sum(!sim$substitutions$synonymous)),
    sprintf("outputs:           alignment_codon.fasta, alignment_amino_acid.fasta,"),
    sprintf("                   substitutions.tsv, branch_modes.tsv%s, run_info.txt",
            ifelse(widget_ok, ", substitution_history.html", ""))),
    out("run_info.txt"))

message(sprintf("done: %d substitutions at %d sites; outputs in %s/",
                nrow(sim$substitutions), n_sites, opts$out))
