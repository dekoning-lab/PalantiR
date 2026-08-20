# MIN7 -- validate before dereferencing, and make has_class() total.
#
# Before the fix: has_class() did `std::string c = a.attr("class");`, which
# THREW on any object without a class attribute instead of answering FALSE; the
# three simulators dereferenced models[0]$type before any validation, so an
# empty model list gave an out-of-bounds error and a plain list gave Rcpp's
# "Object was created without names".
# EXPECTED against the OLD library: FAIL -- the messages are Rcpp internals that
#   name neither the argument nor its expected class, and the empty-list case
#   fails on an index rather than on the list being empty.
# EXPECTED against the FIXED library: PASS -- every malformed argument produces a
#   message naming the argument and what was expected.

suppressMessages(library(PalantiR))

ok <- FALSE
why <- ""

errs <- function(expr, must_mention) {
    m <- tryCatch({ force(expr); NA_character_ },
                  error = function(e) conditionMessage(e))
    if (is.na(m)) return(FALSE)
    all(vapply(must_mention, function(w) grepl(w, m, fixed = TRUE), logical(1)))
}

tryCatch({
    tree_path  <- system.file("extdata", "primates.newick", package = "PalantiR")
    guide_path <- system.file("extdata", "primates_switch.newick", package = "PalantiR")
    tree  <- Phylogeny(tree_path)
    guide <- Phylogeny(guide_path, type = "mode")

    hky <- HasegawaKishinoYano(rep(0.25, 4), transition_rate = 2)
    gtr <- GeneralTimeReversible(c(0.5, 0.5), matrix(c(0, 1, 1, 0), 2, 2))
    psi <- equilibrium_to_fitness(rep(1/20, 20), 1000)
    m   <- MutationSelection(1000, 1, hky, psi, scaling_type = "non-synonymous")
    nuc_seq <- Sequence("ACGTACGT", type = "nucleotide")
    cod_seq <- Sequence("ATGGCC", type = "codon")

    # A classless list, the object that used to make has_class() throw.
    classless <- list(newick = "(a:1,b:1);", type = "phylogeny", n_nodes = 3)

    checks <- c(
        # has_class() is total: a classless list is simply not a Phylogeny.
        "classless list as tree" =
            errs(simulate_over_phylogeny(classless, hky, nuc_seq), c("tree", "Phylogeny")),
        "classless list as model" =
            errs(simulate_over_phylogeny(tree, list(a = 1), nuc_seq),
                 c("substitution_model", "SubstitutionModel")),
        "classless list as sequence" =
            errs(simulate_over_phylogeny(tree, hky, list(a = 1)), c("sequence", "Sequence")),
        # An empty model list must say so, not fail on an index.
        "empty models (interval simulator)" =
            errs(simulate_over_interval_phylogeny(tree, guide, list(), cod_seq, 1),
                 c("substitution_models", "at least one")),
        "empty models (shared substitution)" =
            errs(simulate_with_shared_substitution_heterogeneity(
                     tree, hky, list(), cod_seq, 0), c("substitution_models", "at least one")),
        "empty models (shared time)" =
            errs(simulate_with_shared_time_heterogeneity(
                     tree, hky, list(), cod_seq, 0), c("substitution_models", "at least one")),
        "empty models (Markov modulated)" =
            errs(MarkovModulatedMutationSelection(list(), gtr),
                 c("mutation_selection_models", "at least one")),
        # A list of non-models must be caught before its `type` is read.
        "non-models in the list" =
            errs(simulate_over_interval_phylogeny(tree, guide, list(list(a = 1)), cod_seq, 1),
                 c("substitution_models", "SubstitutionModel"))
    )
    bad <- names(checks)[!checks]

    # Positive control: a well-formed call still gets through the new gauntlet.
    set_palantir_seed(3)
    sim <- simulate_over_phylogeny(tree, hky, nuc_seq, rate = 1)

    if (length(bad) > 0) {
        why <- paste("wrong or missing error for:", paste(bad, collapse = "; "))
    } else if (is.null(sim$alignment) || nrow(sim$alignment) == 0) {
        why <- "positive control: a valid simulate_over_phylogeny() call produced no alignment"
    } else {
        ok <- TRUE
    }
}, error = function(e) why <<- conditionMessage(e))

cat(if (ok) "MIN7: PASS\n" else paste0("MIN7: FAIL ", why, "\n"))
