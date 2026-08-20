# D3: the README demonstration does not run.
#
# Old behaviour: the demo opened with Phylogeny("data/mammals.newick") and
# read.csv("data/amino_acid_fitness_N_1000.csv"). No mammals.newick has ever
# been in data/ (the trees are in inst/extdata), so the demo died on its first
# line with "Could not find file data/mammals.newick"; and the fitness csv has
# now moved to inst/extdata as well. It also plotted sites = 1:10, which -- the
# site index being 0-based -- draws sites two to eleven. This script FAILS
# against that README (error while evaluating the extracted demo).
#
# New behaviour: the demo reads the tree with system.file("extdata",
# "mammals.newick", package = "PalantiR"), takes the fitnesses from the packaged
# amino_acid_fitness_N_1000 dataset, and plots sites = 0:9. Every R chunk of the
# "Short demonstration" section runs to completion. This script PASSES.
#
# The demo is not transcribed here: it is extracted from README.md at run time,
# so the test cannot drift away from the file it is checking.

suppressMessages(library(PalantiR))

# Resolve the repository root by looking for src/Palantir_Core upward from the
# working directory; source-inspection tests SKIP when run without the source
# tree (for example against an installed package).
find_repo <- function() {
    for (cand in c("..", "../..", ".", "../../..")) {
        if (dir.exists(file.path(cand, "src", "Palantir_Core"))) return(normalizePath(cand))
    }
    NA
}

repo <- find_repo()
if (is.na(repo)) { cat("D3: SKIP (source tree not available)\n"); quit(save = "no", status = 0) }

fail <- function(reason) {
    cat("D3: FAIL", reason, "\n")
    quit(save = "no", status = 0)
}

readme <- file.path(repo, "README.md")
if (!file.exists(readme)) fail("README.md not found")
lines <- readLines(readme, warn = FALSE)

# Everything from the "Short demonstration" heading onwards.
start <- grep("^#+\\s*Short demonstration", lines)
if (length(start) != 1L) fail("could not locate the 'Short demonstration' heading")
lines <- lines[start:length(lines)]

# Collect the bodies of the ```R fenced blocks.
fences <- grep("^```", lines)
if (length(fences) %% 2L != 0L) fail("unbalanced code fences in README.md")
code <- character(0)
for (i in seq(1, length(fences), by = 2L)) {
    open <- fences[i]; close <- fences[i + 1L]
    if (!grepl("^```[Rr]\\s*$", lines[open])) next
    if (close > open + 1L) code <- c(code, lines[(open + 1L):(close - 1L)])
}
if (length(code) == 0L) fail("no R code blocks found in the demonstration section")
if (!any(grepl("simulate_over_phylogeny", code)))
    fail("the extracted demonstration does not call simulate_over_phylogeny")

# The demo writes PalantiR_ms.fa into the working directory; run it from the
# repository root, as a reader with a clone would, and clean up afterwards.
old <- setwd(repo)
on.exit({ unlink(file.path(repo, "PalantiR_ms.fa")); setwd(old) }, add = TRUE)

env <- new.env(parent = globalenv())
res <- tryCatch({
    eval(parse(text = paste(code, collapse = "\n")), envir = env)
    TRUE
}, error = function(e) conditionMessage(e))
if (!isTRUE(res)) fail(paste("the README demonstration failed:", res))

if (!exists("sim", envir = env, inherits = FALSE)) fail("the demo produced no `sim`")
sim <- get("sim", envir = env)
if (!inherits(sim, "Simulation")) fail("`sim` is not a Simulation")
if (ncol(sim$alignment) != 100L)
    fail(paste("the simulated alignment has", ncol(sim$alignment), "sites, expected 100"))
if (!file.exists(file.path(repo, "PalantiR_ms.fa")))
    fail("the demo's as.fasta call wrote no file")

# The plot call must ask for the first ten sites, which are 0:9.
if (!any(grepl("sites\\s*=\\s*0:9", code)))
    fail("the demo still plots 1-based sites; the site index is 0-based")

cat("D3: PASS\n")
