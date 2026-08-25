#!/usr/bin/env Rscript
#
# History-level validation of branch-heterogeneous Goldman--Yang 1994 models.
#
# This diagnostic simulates a deep 50-taxon tree under synonymous scaling and
# estimates dN/dS separately in its two principal subtrees from the complete
# substitution history.  The estimator uses the Goldman--Yang definition of
# synonymous and nonsynonymous sites, rather than the raw N/S event ratio.
#
# Usage:
#   Rscript validation/gy94_deep_tree_history_dnds.R \
#     --out=runs/gy94-deep-tree-dnds-20260825

suppressMessages(library(PalantiR))
use_genetic_code("Standard nuclear")

defaults <- list(
    out = file.path("runs", "gy94-deep-tree-dnds-20260825"),
    sites = 5000L,
    seed = 20260825L)

args <- commandArgs(trailingOnly = TRUE)
if(any(args %in% c("-h", "--help"))) {
    cat(paste0(
        "Options:\n",
        "  --out=DIR   artifact directory (", defaults$out, ")\n",
        "  --sites=N   codon sites (", defaults$sites, ")\n",
        "  --seed=N    PalantiR RNG seed (", defaults$seed, ")\n"))
    quit(save = "no", status = 0)
}

opts <- defaults
for(arg in args) {
    if(!grepl("^--(out|sites|seed)=", arg)) stop("Unknown argument: ", arg)
    key <- sub("=.*$", "", sub("^--", "", arg))
    opts[[key]] <- sub("^[^=]*=", "", arg)
}
opts$sites <- as.integer(opts$sites)
opts$seed <- as.integer(opts$seed)
if(is.na(opts$sites) || opts$sites < 1L) stop("--sites should be positive")
if(is.na(opts$seed)) stop("--seed should be an integer")

if(dir.exists(opts$out) && length(list.files(opts$out, all.files = TRUE,
                                             no.. = TRUE))) {
    stop("Refusing to overwrite non-empty artifact directory: ", opts$out)
}
dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
out <- function(filename) file.path(opts$out, filename)

# ---------------------------------------------------------------------------
# Parameters and a deterministic, ultrametric 50-tip topology
# ---------------------------------------------------------------------------

n_sites <- opts$sites
seed <- opts$seed
kappa <- 4
omega <- c(`0` = 0.1, `1` = 1.0)
nucleotide_frequencies <- c(T = 0.22, C = 0.28, A = 0.27, G = 0.23)

# In synonymous scaling, branch length is measured in expected synonymous
# changes per codon.  Every root-to-tip path is therefore exactly 20 in that
# currency.  Most depth is placed on the two basal stems, keeping the complete
# event-history artifact tractable while retaining a fully bifurcating crown.
root_to_tip_depth <- 20
basal_stem <- 19
crown_depth <- root_to_tip_depth - basal_stem
internal_edge <- 0.15

render_clade <- function(labels, mode, depth = 0) {
    if(length(labels) < 2L) stop("render_clade requires at least two labels")
    split_at <- floor(length(labels) / 2L)
    groups <- list(labels[seq_len(split_at)], labels[-seq_len(split_at)])

    render_child <- function(group) {
        if(length(group) == 1L) {
            terminal <- crown_depth - depth
            if(terminal <= 0) stop("Crown construction exhausted its depth")
            return(list(
                tree = sprintf("%s:%.15g", group, terminal),
                modes = sprintf("%s:%d", group, mode)))
        }
        child <- render_clade(group, mode, depth + internal_edge)
        list(
            tree = sprintf("%s:%.15g", child$tree, internal_edge),
            modes = sprintf("%s:%d", child$modes, mode))
    }

    children <- lapply(groups, render_child)
    list(
        tree = sprintf("(%s,%s)", children[[1]]$tree, children[[2]]$tree),
        modes = sprintf("(%s,%s)", children[[1]]$modes, children[[2]]$modes))
}

left <- render_clade(sprintf("L%02d", seq_len(25L)), mode = 0L)
right <- render_clade(sprintf("R%02d", seq_len(25L)), mode = 1L)
tree_newick <- sprintf("(%s:%.15g,%s:%.15g);",
                       left$tree, basal_stem, right$tree, basal_stem)
mode_newick <- sprintf("(%s:0,%s:1);", left$modes, right$modes)
writeLines(tree_newick, out("tree_50_taxa_depth20.newick"))
writeLines(mode_newick, out("mode_tree_left_omega0.1_right_omega1.newick"))

tree <- Phylogeny(out("tree_50_taxa_depth20.newick"))
mode_tree <- Phylogeny(out("mode_tree_left_omega0.1_right_omega1.newick"),
                       type = "mode")
branch_table <- phylogeny_to_intervals(tree, mode_tree)
branch_table$branch_length <- branch_table$to - branch_table$from
branch_table$subtree <- ifelse(branch_table$state == 0,
                               "left_omega_0.1", "right_omega_1.0")
if(nrow(branch_table) != 98L || length(unique(branch_table$node)) != 98L) {
    stop("The generated 50-tip binary tree should contain exactly 98 branches")
}
write.table(branch_table, out("branch_modes_and_lengths.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# GY94 model and neutral synonymous/nonsynonymous opportunities
# ---------------------------------------------------------------------------

pi <- F1x4(nucleotide_frequencies)
branch_process <- GY94BranchModel(
    models = list(
        GY94(omega = omega[["0"]], kappa = kappa, frequencies = pi,
             scaling_type = "synonymous"),
        GY94(omega = omega[["1"]], kappa = kappa, frequencies = pi,
             scaling_type = "synonymous")),
    mode_phylogeny = mode_tree,
    start_mode = 0L)
site_model <- GY94SiteModel(list(branch_process = branch_process),
                            n_sites = n_sites)

genetic_code <- GeneticCode()
codons <- names(genetic_code)
amino_acids <- unname(genetic_code)
codon_chars <- strsplit(codons, "")
one_step <- outer(seq_along(codons), seq_along(codons),
                  Vectorize(function(i, j) {
                      i != j && sum(codon_chars[[i]] != codon_chars[[j]]) == 1L
                  }))
same_amino_acid <- outer(amino_acids, amino_acids, "==")
synonymous_mask <- one_step & same_amino_acid
nonsynonymous_mask <- one_step & !same_amino_acid

# At omega=1, stationary synonymous and nonsynonymous fluxes determine the
# Goldman--Yang counts of synonymous and nonsynonymous sites per codon:
# S = 3 rS/(rS+rN), N = 3 rN/(rS+rN).  Standard scaling is convenient here
# because rS+rN=1; the opportunity proportions are invariant to common scale.
neutral <- GY94(omega = 1, kappa = kappa, frequencies = pi,
                scaling_type = "standard")
stationary_flux <- function(model, mask) {
    sum(as.numeric(model$equilibrium) * rowSums(model$transition * mask))
}
rS0 <- stationary_flux(neutral, synonymous_mask)
rN0 <- stationary_flux(neutral, nonsynonymous_mask)
pS <- rS0 / (rS0 + rN0)
pN <- rN0 / (rS0 + rN0)
S_sites <- 3 * pS
N_sites <- 3 * pN

opportunities <- data.frame(
    event_class = c("synonymous", "nonsynonymous"),
    neutral_stationary_flux = c(rS0, rN0),
    neutral_opportunity_fraction = c(pS, pN),
    gy94_sites_per_codon = c(S_sites, N_sites))
write.table(opportunities, out("goldman_yang_opportunities.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Confirm the matrices actually passed to the simulator have the expected
# synonymous-scaled rates.  This uses the post-mixture-normalization models.
effective_models <- site_model$models[[1]]$models
matrix_rates <- do.call(rbind, lapply(seq_along(effective_models), function(i) {
    model <- effective_models[[i]]
    data.frame(
        mode = i - 1L,
        assigned_omega = omega[[as.character(i - 1L)]],
        synonymous_rate = stationary_flux(model, synonymous_mask),
        nonsynonymous_rate = stationary_flux(model, nonsynonymous_mask),
        total_rate = stationary_flux(model, synonymous_mask | nonsynonymous_mask),
        scaling_denominator = model$scaling)
}))
write.table(matrix_rates, out("effective_matrix_rates.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
if(any(abs(matrix_rates$synonymous_rate - 1) > 1e-11)) {
    stop("A synonymous-scaled GY94 matrix did not have synonymous rate one")
}

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

set.seed(seed)
set_palantir_seed(seed)
message(sprintf("Simulating %d codon sites on a 50-tip depth-%g tree ...",
                n_sites, root_to_tip_depth))
started <- Sys.time()
simulation <- simulate_gy94_site_model(tree, site_model,
                                       rescale_method = "exact")
elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

if(nrow(simulation$alignment) != 50L || ncol(simulation$alignment) != n_sites) {
    stop("The simulated alignment does not have the requested dimensions")
}
as.fasta(simulation$alignment, file = out("alignment_codon.fasta"))

history <- simulation$substitutions
mode_by_node <- setNames(branch_table$state, as.character(branch_table$node))
history$mode <- unname(mode_by_node[as.character(history$node)])
if(anyNA(history$mode)) stop("At least one history event could not be mapped to a branch mode")
history$subtree <- ifelse(history$mode == 0,
                          "left_omega_0.1", "right_omega_1.0")

history_connection <- gzfile(out("substitution_history.tsv.gz"), open = "wt",
                             compression = 6)
write.table(history, history_connection, sep = "\t", quote = FALSE,
            row.names = FALSE)
close(history_connection)

# ---------------------------------------------------------------------------
# History-level Goldman--Yang dN/dS
# ---------------------------------------------------------------------------

summaries <- vector("list", length(omega))
per_site <- vector("list", length(omega))
for(mode in as.integer(names(omega))) {
    label <- if(mode == 0L) "left_omega_0.1" else "right_omega_1.0"
    events <- history$mode == mode
    synonymous_by_site <- tabulate(history$site[events & history$synonymous] + 1L,
                                   nbins = n_sites)
    nonsynonymous_by_site <- tabulate(history$site[events & !history$synonymous] + 1L,
                                      nbins = n_sites)
    per_site[[mode + 1L]] <- data.frame(
        site = seq_len(n_sites) - 1L,
        mode = mode,
        subtree = label,
        synonymous_events = synonymous_by_site,
        nonsynonymous_events = nonsynonymous_by_site)

    C_S <- sum(synonymous_by_site)
    C_N <- sum(nonsynonymous_by_site)
    dS_by_site <- synonymous_by_site / S_sites
    dN_by_site <- nonsynonymous_by_site / N_sites
    empirical_ratio <- mean(dN_by_site) / mean(dS_by_site)

    # Sitewise delta-method uncertainty.  Sites are independent simulation
    # replicates; substitutions within a site and across branches need not be
    # treated as independent Poisson events.
    influence <- (dN_by_site - empirical_ratio * dS_by_site) / mean(dS_by_site)
    ratio_se <- stats::sd(influence) / sqrt(n_sites)
    log_se <- ratio_se / empirical_ratio
    ratio_ci <- empirical_ratio * exp(c(-1, 1) * stats::qnorm(0.975) * log_se)

    branches <- branch_table$state == mode
    total_branch_length <- sum(branch_table$branch_length[branches])
    exposure <- n_sites * total_branch_length
    model_rates <- matrix_rates[matrix_rates$mode == mode, ]

    summaries[[mode + 1L]] <- data.frame(
        mode = mode,
        subtree = label,
        assigned_omega = omega[[as.character(mode)]],
        tips = 25L,
        branches = sum(branches),
        root_to_tip_synonymous_depth = root_to_tip_depth,
        total_branch_length = total_branch_length,
        codon_branch_exposure = exposure,
        synonymous_events = C_S,
        nonsynonymous_events = C_N,
        raw_n_over_s = C_N / C_S,
        expected_raw_n_over_s = omega[[as.character(mode)]] * rN0 / rS0,
        empirical_synonymous_rate = C_S / exposure,
        expected_synonymous_rate = model_rates$synonymous_rate,
        empirical_nonsynonymous_rate = C_N / exposure,
        expected_nonsynonymous_rate = model_rates$nonsynonymous_rate,
        gy94_dS_tree_length = C_S / (n_sites * S_sites),
        gy94_dN_tree_length = C_N / (n_sites * N_sites),
        empirical_dN_dS = empirical_ratio,
        dN_dS_sitewise_se = ratio_se,
        dN_dS_95pct_lower = ratio_ci[[1]],
        dN_dS_95pct_upper = ratio_ci[[2]],
        relative_error = empirical_ratio / omega[[as.character(mode)]] - 1)
}
summary_table <- do.call(rbind, summaries)
per_site_table <- do.call(rbind, per_site)
write.table(summary_table, out("history_dnds_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
per_site_connection <- gzfile(out("history_counts_by_site.tsv.gz"), open = "wt",
                              compression = 6)
write.table(per_site_table, per_site_connection, sep = "\t", quote = FALSE,
            row.names = FALSE)
close(per_site_connection)

# Preserve compact R objects needed for exact downstream recalculation without
# duplicating the large, language-neutral substitution-history table.
saveRDS(list(
    parameters = list(n_sites = n_sites, seed = seed, kappa = kappa,
                      omega = omega,
                      nucleotide_frequencies = nucleotide_frequencies,
                      scaling_type = "synonymous",
                      root_to_tip_synonymous_depth = root_to_tip_depth),
    codon_frequencies = pi,
    site_model = site_model,
    branch_table = branch_table,
    opportunities = opportunities,
    matrix_rates = matrix_rates,
    summary = summary_table,
    alignment = simulation$alignment),
    out("diagnostic_objects.rds"), compress = "xz")

git_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]],
    error = function(e) NA_character_)
writeLines(c(
    "PalantiR GY94 deep-tree history-level dN/dS diagnostic",
    sprintf("completed: %s", format(Sys.time(), tz = "America/Edmonton")),
    sprintf("elapsed_seconds: %.3f", elapsed_seconds),
    sprintf("PalantiR_version: %s", as.character(packageVersion("PalantiR"))),
    sprintf("PalantiR_git_commit: %s", git_commit),
    sprintf("R_version: %s", R.version.string),
    sprintf("genetic_code: %s", get_genetic_code()),
    sprintf("taxa: %d", nrow(simulation$alignment)),
    sprintf("codon_sites: %d", n_sites),
    sprintf("seed: %d", seed),
    sprintf("kappa: %.10g", kappa),
    sprintf("F1x4_nucleotide_frequencies_T_C_A_G: %s",
            paste(format(nucleotide_frequencies, digits = 10), collapse = ",")),
    "scaling: synonymous",
    sprintf("root_to_tip_depth_in_synonymous_scaled_branch_units: %.10g",
            root_to_tip_depth),
    sprintf("neutral_synonymous_opportunity_fraction: %.12g", pS),
    sprintf("neutral_nonsynonymous_opportunity_fraction: %.12g", pN),
    sprintf("Goldman_Yang_synonymous_sites_per_codon: %.12g", S_sites),
    sprintf("Goldman_Yang_nonsynonymous_sites_per_codon: %.12g", N_sites),
    "",
    "Goldman--Yang history estimator:",
    "  S_sites = 3 rS0/(rS0+rN0); N_sites = 3 rN0/(rS0+rN0),",
    "  where rS0 and rN0 are stationary synonymous and nonsynonymous fluxes",
    "  from the omega=1 generator with the fitted kappa and codon frequencies.",
    "  dS = C_S/(codon_sites*S_sites); dN = C_N/(codon_sites*N_sites).",
    "  Reported dN/dS is dN/dS, not the raw nonsynonymous/synonymous count ratio.",
    "  This is the construction in PAML codeml eigenQcodon(), codeml.c lines",
    "  3323-3366: https://github.com/abacus-gene/paml/blob/master/src/codeml.c",
    "",
    capture.output(print(summary_table, row.names = FALSE, digits = 8))),
    out("run_info.txt"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if(length(script_arg)) {
    script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
    invisible(file.copy(script_path, out("gy94_deep_tree_history_dnds.R"),
                        overwrite = TRUE))
}

artifact_files <- list.files(opts$out, full.names = TRUE)
manifest <- data.frame(
    file = basename(artifact_files),
    bytes = unname(file.info(artifact_files)$size),
    md5 = unname(tools::md5sum(artifact_files)),
    stringsAsFactors = FALSE)
write.table(manifest, out("MANIFEST.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

cat("\nCompleted in", sprintf("%.1f", elapsed_seconds), "seconds.\n")
print(summary_table[, c("subtree", "assigned_omega", "synonymous_events",
                        "nonsynonymous_events", "raw_n_over_s",
                        "empirical_dN_dS", "dN_dS_sitewise_se",
                        "dN_dS_95pct_lower", "dN_dS_95pct_upper",
                        "relative_error")],
      row.names = FALSE, digits = 8)
cat("\nArtifact:", normalizePath(opts$out), "\n")
