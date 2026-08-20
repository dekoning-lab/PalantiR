# M10 -- join() offsets the 0-based `site` column and normalizes `intervals`.
#
# The `site` column of a substitution table is a 0-based index into that
# simulation's OWN alignment. join() rbind()ed the tables untouched while
# cbind()ing the alignments, so site 0 of the second simulation landed on top of
# site 0 of the first: a joined simulation of two 2-site runs claimed 4 alignment
# columns but only sites 0 and 1, with unrelated substitutions overlaid. It also
# concatenated the `intervals` members verbatim, one entry per SIMULATION where
# every other member is per SITE.
#
# EXPECTED against the OLD library: FAIL -- the joined substitution table has
#   sites 0,1 for a 4-column alignment, and joining two interval simulations
#   gives an unnamed 2-element interval list.
# EXPECTED against the FIXED library: PASS -- sites are 0,1,2,3, the joined table
#   is exactly the two tables with the second shifted by ncol(first alignment),
#   the alignment columns line up with those indices, and the joined intervals
#   are one table per site named by the site index.

suppressMessages(library(PalantiR))

ID <- "M10"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

bare <- function(df) { rownames(df) <- NULL; df }

result <- try({
    mk <- function(txt, type = "phylogeny") {
        f <- tempfile(fileext = ".newick"); writeLines(txt, f); Phylogeny(f, type = type)
    }
    tree <- mk("(((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0);")

    hky <- HasegawaKishinoYano(equilibrium = c(.25, .25, .25, .25))
    set.seed(23)
    m1 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))
    m2 <- MutationSelection(1000, 1e-8, hky, 1 + rnorm(20, 0, 5e-4))

    a <- simulate_over_phylogeny(tree, m1, sample_sequence(model = m1, length = 2), rate = 1)
    b <- simulate_over_phylogeny(tree, m2, sample_sequence(model = m2, length = 2), rate = 1)
    if (nrow(a$substitutions) == 0 || nrow(b$substitutions) == 0) {
        fail("inconclusive: a simulation produced no substitutions")
    }

    joined <- join(a, b)

    # 1. site indices span the joined alignment
    sites <- sort(unique(joined$substitutions$site))
    if (ncol(joined$alignment) != 4) {
        fail(sprintf("joined alignment has %d columns, expected 4",
                     ncol(joined$alignment)))
    }
    if (!identical(as.numeric(sites), as.numeric(0:3))) {
        fail(sprintf("joined sites are %s, expected 0,1,2,3",
                     paste(sites, collapse = ",")))
    }

    # 2. the joined table IS the two tables, the second shifted by ncol(a)
    shifted <- b$substitutions
    shifted$site <- shifted$site + ncol(a$alignment)
    expected <- bare(rbind(a$substitutions, shifted))
    if (!isTRUE(all.equal(bare(joined$substitutions), expected))) {
        fail("the joined substitution table is not the two tables with the second offset")
    }

    # 3. those indices address the alignment columns they came from
    A <- unclass(a$alignment); B <- unclass(b$alignment); J <- unclass(joined$alignment)
    if (!all(J[, 1:2] == A) || !all(J[, 3:4] == B)) {
        fail("the joined alignment columns are not the two alignments in order")
    }
    if (!identical(rownames(J), rownames(A))) {
        fail("the joined alignment lost its taxon names")
    }

    # 4. homogeneous simulations contribute no intervals at all
    if (!is.null(joined$intervals)) {
        fail(sprintf("joining two homogeneous simulations produced intervals (%s)",
                     paste(class(joined$intervals), collapse = "/")))
    }

    # 5. interval simulations give one table per site, named by the site index
    mode_tree <- mk("(((A:0,B:0):0,(C:1,D:1):1):0);", type = "mode")
    ia <- simulate_over_interval_phylogeny(
        tree, mode_tree, list(m1, m2), sample_sequence(model = m1, length = 2),
        start_mode = 0)
    ib <- simulate_over_interval_phylogeny(
        tree, mode_tree, list(m1, m2), sample_sequence(model = m2, length = 2),
        start_mode = 0)
    ij <- join(ia, ib)

    if (!is.list(ij$intervals) || is.data.frame(ij$intervals)) {
        fail("joined intervals are not a list")
    }
    if (length(ij$intervals) != ncol(ij$alignment)) {
        fail(sprintf("joined intervals have %d entries for %d alignment columns",
                     length(ij$intervals), ncol(ij$alignment)))
    }
    if (!identical(names(ij$intervals), as.character(0:3))) {
        fail(sprintf("joined intervals are named %s, expected the site indices 0..3",
                     paste(names(ij$intervals), collapse = ",")))
    }
    if (!all(vapply(ij$intervals, is.data.frame, logical(1)))) {
        fail("a joined interval entry is not an interval table")
    }
    if (!isTRUE(all.equal(bare(ij$intervals[["0"]]), bare(ia$intervals))) ||
        !isTRUE(all.equal(bare(ij$intervals[["3"]]), bare(ib$intervals)))) {
        fail("joined interval tables are not the source simulations' tables")
    }
    if (!identical(as.numeric(sort(unique(ij$substitutions$site))), as.numeric(0:3))) {
        fail("interval simulations were joined without offsetting their sites")
    }

    pass()
}, silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
