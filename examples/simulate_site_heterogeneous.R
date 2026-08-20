#!/usr/bin/env Rscript
#
# simulate_site_heterogeneous.R
#
# Command-line driver for a time-homogeneous mutation-selection simulation in
# which sites differ in their amino-acid fitness profiles. Each site is
# assigned one of a set of 20-dimensional amino-acid equilibrium profiles;
# each profile is converted to a fitness vector at the chosen population size
# with equilibrium_to_fitness(); sites sharing a profile are simulated
# together under one MutationSelection model and the per-profile simulations
# are combined with join().
#
# Usage:
#   Rscript simulate_site_heterogeneous.R [--option=value ...]
#
# Options (defaults in parentheses):
#   --sites=N             number of alignment sites to simulate (100)
#   --population-size=N   effective population size (5000)
#   --mutation-rate=X     per-generation mutation rate (1e-8)
#   --tree=FILE           newick tree with branch lengths
#                         (mammals.newick shipped with the package)
#   --profiles=FILE       whitespace-delimited table of amino-acid equilibrium
#                         frequencies, one profile per row, 20 columns in the
#                         order A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
#                         (the 50 profiles in the packaged psi_c50_1 data set)
#   --seed=N              random seed (20260820)
#   --out=DIR             output directory (site_heterogeneous_output)
#
# Outputs, written to --out:
#   alignment_codon.fasta       simulated alignment, codon sequences
#   alignment_amino_acid.fasta  the same alignment, translated
#   substitutions.tsv           full substitution history (one row per event)
#   site_profiles.tsv           profile index assigned to each site
#   substitution_history.html   interactive substitution-history plot
#   run_info.txt                parameters, seed, versions, output summary
#
# Example:
#   Rscript simulate_site_heterogeneous.R --sites=200 --population-size=10000 --out=run1

suppressMessages(library(PalantiR))

# ---------------------------------------------------------------------------
# Command-line options
# ---------------------------------------------------------------------------

defaults <- list(
    sites           = "100",
    population_size = "5000",
    mutation_rate   = "1e-8",
    tree            = system.file("extdata", "mammals.newick", package = "PalantiR"),
    profiles        = "",          # empty: use the packaged psi_c50_1 data set
    seed            = "20260820",
    out             = "site_heterogeneous_output")

argv <- commandArgs(trailingOnly = TRUE)
if (any(argv %in% c("-h", "--help"))) {
    cat("Options (--option=value):\n")
    for (k in names(defaults))
        cat(sprintf("  --%-16s (default: %s)\n", gsub("_", "-", k),
                    ifelse(defaults[[k]] == "", "packaged psi_c50_1", defaults[[k]])))
    quit(status = 0)
}
opts <- defaults
for (a in argv) {
    if (!grepl("^--[a-z-]+=.+", a)) stop("unrecognized argument: ", a, " (see --help)")
    key <- gsub("-", "_", sub("=.*", "", sub("^--", "", a)))
    if (!key %in% names(opts)) stop("unknown option in: ", a, " (see --help)")
    opts[[key]] <- sub("^[^=]*=", "", a)
}

n_sites  <- as.integer(opts$sites)
pop_size <- as.numeric(opts$population_size)
mu       <- as.numeric(opts$mutation_rate)
seed     <- as.integer(opts$seed)
stopifnot(n_sites > 0, pop_size > 0, mu > 0)

# Seed both random number generators: PalantiR's C++ generator drives the
# simulation itself, R's generator draws the site-to-profile assignment.
set_palantir_seed(seed)
set.seed(seed)

# ---------------------------------------------------------------------------
# Inputs: tree and fitness profiles
# ---------------------------------------------------------------------------

phylogeny <- Phylogeny(opts$tree)

if (nzchar(opts$profiles)) {
    profiles <- as.matrix(read.table(opts$profiles))
    if (ncol(profiles) != 20)
        stop("--profiles must have 20 columns (amino-acid frequencies); got ",
             ncol(profiles))
} else {
    data(psi_c50_1, package = "PalantiR", envir = environment())
    profiles <- as.matrix(psi_c50_1)
}
n_profiles <- nrow(profiles)

# ---------------------------------------------------------------------------
# Site-to-profile assignment and per-profile models
# ---------------------------------------------------------------------------

# Sites draw their profile with replacement. Sites that share a profile are
# simulated together as one block, so in the output alignment they appear
# consecutively; site_profiles.tsv records the assignment site by site.
assignment <- sort(sample.int(n_profiles, n_sites, replace = TRUE))
used <- unique(assignment)

nucleotide_model <- HasegawaKishinoYano(equilibrium = rep(0.25, 4))

models <- lapply(used, function(p) {
    psi <- equilibrium_to_fitness(profiles[p, ], pop_size)
    MutationSelection(
        population_size  = pop_size,
        mutation_rate    = mu,
        nucleotide_model = nucleotide_model,
        fitness          = psi)
})

# ---------------------------------------------------------------------------
# Simulation: one block per profile, combined with join()
# ---------------------------------------------------------------------------

message(sprintf("simulating %d sites over %d fitness profiles (N = %g) ...",
                n_sites, length(used), pop_size))

block_sims <- lapply(seq_along(used), function(i) {
    n_block <- sum(assignment == used[i])
    root    <- sample_sequence(models[[i]], n_block)
    simulate_over_phylogeny(
        phylogeny = phylogeny,
        model     = models[[i]],
        sequence  = root)
})
sim <- if (length(block_sims) == 1) block_sims[[1]] else do.call(join, block_sims)

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

write.table(data.frame(site = seq_len(n_sites) - 1L, profile = assignment),
            out("site_profiles.tsv"),
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
    "PalantiR site-heterogeneous simulation",
    sprintf("date:              %s", format(Sys.time())),
    sprintf("PalantiR version:  %s", as.character(packageVersion("PalantiR"))),
    sprintf("R version:         %s", R.version.string),
    sprintf("tree:              %s", opts$tree),
    sprintf("profiles:          %s (%d profiles, %d used)",
            ifelse(nzchar(opts$profiles), opts$profiles, "packaged psi_c50_1"),
            n_profiles, length(used)),
    sprintf("sites:             %d", n_sites),
    sprintf("population size:   %g", pop_size),
    sprintf("mutation rate:     %g", mu),
    sprintf("seed:              %d", seed),
    sprintf("substitutions:     %d events (%d non-synonymous)",
            nrow(sim$substitutions), sum(!sim$substitutions$synonymous)),
    sprintf("outputs:           alignment_codon.fasta, alignment_amino_acid.fasta,"),
    sprintf("                   substitutions.tsv, site_profiles.tsv%s, run_info.txt",
            ifelse(widget_ok, ", substitution_history.html", ""))),
    out("run_info.txt"))

message(sprintf("done: %d substitutions at %d sites; outputs in %s/",
                nrow(sim$substitutions), n_sites, opts$out))
