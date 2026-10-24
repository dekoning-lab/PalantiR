# D5: the C50 profiles have V and W in the wrong columns.
#
# Old behaviour: psi_c50_1 (Table 4 of the supplement to Le, Gascuel &
# Lartillot 2008) held the table's W frequencies in column 18, documented as V,
# and its V frequencies in column 19, documented as W. Profile 2, 41% I and 43%
# V in the paper, read as 41% I and 43% W. This script FAILS against that.
#
# New behaviour (1.3.0): psi_c50_1 and psi_c50_2 (Table 5) carry amino-acid
# column names in the package's order, and multiplying a row by its printed
# total recovers the table's integers. Spot checks below cover the two
# aromatic-rich profiles and the isoleucine/valine profile, and two composition
# checks guard against any future V/W exchange across all 50 rows. The IQ-TREE
# version, psi_c50_iqtree, is pinned to IQ-TREE's C50pi1, which is listed in
# ARNDCQEGHILKMFPSTWYV order and so checks the reordering. This script PASSES.

suppressMessages(library(PalantiR))

fail <- function(reason) {
    cat("D5: FAIL", reason, "\n")
    quit(save = "no", status = 0)
}

aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

env <- new.env()
data(list = c("psi_c50_1", "psi_c50_2", "psi_c50_iqtree"), package = "PalantiR", envir = env)

# Printed integers (whole percentages) and row totals from Tables 4 and 5.
spots <- list(
    list(set = "psi_c50_1", profile = 2,  total = 99, pct = c(F = 1,  I = 41, V = 43, W = 0, Y = 0)),
    list(set = "psi_c50_1", profile = 5,  total = 98, pct = c(F = 41, I = 1,  V = 1,  W = 4, Y = 44)),
    list(set = "psi_c50_1", profile = 39, total = 99, pct = c(F = 22, I = 4,  V = 5,  W = 9, Y = 22)),
    list(set = "psi_c50_2", profile = 2,  total = 99, pct = c(F = 1,  I = 40, V = 43, W = 0, Y = 0)),
    list(set = "psi_c50_2", profile = 5,  total = 98, pct = c(F = 40, I = 1,  V = 1,  W = 5, Y = 43))
)

for (s in spots) {
    ps <- get(s$set, envir = env)
    if (!identical(colnames(ps), aa))
        fail(paste(s$set, "does not have the 20 amino-acid column names"))
    got <- ps[s$profile, names(s$pct)] * s$total
    if (max(abs(got - s$pct)) > 1e-6)
        fail(sprintf("%s profile %d: %s, expected %s from the table", s$set, s$profile,
                     paste(names(got), round(got, 2), collapse = " "),
                     paste(names(s$pct), s$pct, collapse = " ")))
}

# Composition: V is common and tracks I; W is rare. Exchanging the two columns
# breaks both.
for (nm in c("psi_c50_1", "psi_c50_2")) {
    ps <- get(nm, envir = env)
    if (mean(ps[, "V"]) < 5 * mean(ps[, "W"]))
        fail(sprintf("%s: mean V %.3f is not well above mean W %.3f", nm,
                     mean(ps[, "V"]), mean(ps[, "W"])))
    if (cor(ps[, "V"], ps[, "I"]) < 0.5)
        fail(sprintf("%s: V and I are uncorrelated across profiles (r = %.2f)", nm,
                     cor(ps[, "V"], ps[, "I"])))
}

# IQ-TREE's C50pi1 = 0.1357566757 (A) ... 0.2060115780 (S) 0.1011347978 (T)
# 0.0012443033 (W) 0.0056104605 (Y) 0.0093801079 (V), weight 0.0164297003.
iq <- get("psi_c50_iqtree", envir = env)
iq_w <- get("psi_c50_iqtree_weights", envir = env)
if (!is.numeric(iq) || !identical(dim(iq), c(50L, 20L)) || !identical(colnames(iq), aa))
    fail("psi_c50_iqtree is not a 50 x 20 matrix with the amino-acid column names")
iq1 <- c(A = 0.1357566757, S = 0.2060115780, T = 0.1011347978,
         W = 0.0012443033, Y = 0.0056104605, V = 0.0093801079)
if (max(abs(iq[1, names(iq1)] - iq1)) > 1e-12)
    fail("psi_c50_iqtree profile 1 does not match IQ-TREE's C50pi1")
if (max(abs(rowSums(iq) - 1)) > 1e-9)
    fail("psi_c50_iqtree rows do not sum to 1")
if (length(iq_w) != 50 || abs(sum(iq_w) - 1) > 1e-9 || abs(iq_w[1] - 0.0164297003) > 1e-12)
    fail("psi_c50_iqtree_weights are not IQ-TREE's 50 C50 weights")

cat("D5: PASS\n")
