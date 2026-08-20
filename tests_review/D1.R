# D1: data() datasets are unusable.
#
# Old behaviour: data/ held amino_acid_equilibrium.csv,
# amino_acid_fitness_N_1000.csv and psi-c50-1.txt. R's data() reads a bare .csv
# in data/ with read.table(sep = ";"), so data(amino_acid_equilibrium) produced a
# 6 x 1 data frame whose single column held the whole comma-separated line, and
# psi-c50-1 was not a syntactic name at all. This script FAILS against that
# (wrong dimensions / object not found).
#
# New behaviour: data/ holds .rda files. data(amino_acid_equilibrium) and
# data(amino_acid_fitness_N_1000) each yield a 6 x 20 numeric matrix with the 20
# amino-acid column names, and data(psi_c50_1) a 50 x 20 numeric matrix; the raw
# files are in inst/extdata. This script PASSES.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

fail <- function(reason) {
    cat("D1: FAIL", reason, "\n")
    quit(save = "no", status = 0)
}

aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

env <- new.env()
ok <- tryCatch({
    data(list = c("amino_acid_equilibrium", "amino_acid_fitness_N_1000", "psi_c50_1"),
         package = "PalantiR", envir = env)
    TRUE
}, error = function(e) conditionMessage(e))
if (!isTRUE(ok)) fail(paste("data() errored:", ok))

for (nm in c("amino_acid_equilibrium", "amino_acid_fitness_N_1000", "psi_c50_1")) {
    if (!exists(nm, envir = env, inherits = FALSE)) fail(paste("dataset", nm, "did not load"))
}

eq <- get("amino_acid_equilibrium", envir = env)
fi <- get("amino_acid_fitness_N_1000", envir = env)
ps <- get("psi_c50_1", envir = env)

if (!is.numeric(eq) || !identical(dim(eq), c(6L, 20L)))
    fail(paste("amino_acid_equilibrium is", class(eq)[1], "of dim",
               paste(dim(eq), collapse = "x"), "- expected numeric 6x20"))
if (!identical(colnames(eq), aa))
    fail("amino_acid_equilibrium column names are not the 20 amino acids")
if (max(abs(rowSums(eq) - 1)) > 1e-8)
    fail("amino_acid_equilibrium rows do not sum to 1")

if (!is.numeric(fi) || !identical(dim(fi), c(6L, 20L)))
    fail(paste("amino_acid_fitness_N_1000 is", class(fi)[1], "of dim",
               paste(dim(fi), collapse = "x"), "- expected numeric 6x20"))
if (!identical(colnames(fi), aa))
    fail("amino_acid_fitness_N_1000 column names are not the 20 amino acids")
if (max(fi) != 1 || min(fi) <= 0)
    fail("amino_acid_fitness_N_1000 is not a set of relative fitnesses in (0, 1]")

if (!is.numeric(ps) || !identical(dim(ps), c(50L, 20L)))
    fail(paste("psi_c50_1 is", class(ps)[1], "of dim",
               paste(dim(ps), collapse = "x"), "- expected numeric 50x20"))
if (max(abs(rowSums(ps) - 1)) > 1e-8)
    fail("psi_c50_1 rows do not sum to 1")

# The raw files moved to inst/extdata and are still shipped.
for (f in c("amino_acid_equilibrium.csv", "amino_acid_fitness_N_1000.csv", "psi-c50-1.txt")) {
    p <- system.file("extdata", f, package = "PalantiR")
    if (!nzchar(p)) fail(paste("raw file", f, "is not in inst/extdata"))
}

# The packaged matrices agree with the raw files they were built from.
raw_eq <- as.matrix(read.csv(system.file("extdata", "amino_acid_equilibrium.csv",
                                         package = "PalantiR"), check.names = FALSE))
if (!isTRUE(all.equal(eq, raw_eq))) fail("amino_acid_equilibrium differs from its raw csv")
raw_ps <- as.matrix(read.delim(system.file("extdata", "psi-c50-1.txt", package = "PalantiR"),
                               header = FALSE))[, 1:20]
dimnames(raw_ps) <- NULL
if (!isTRUE(all.equal(ps, raw_ps))) fail("psi_c50_1 differs from its raw txt")

cat("D1: PASS\n")
