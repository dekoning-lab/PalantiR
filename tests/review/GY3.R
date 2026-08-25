# GY3 -- arbitrary discrete site assignments, mixed homogeneous/branch classes,
# and propagation of site-class metadata through simulation, joining and widgets.

suppressMessages(library(PalantiR))

ID <- "GY3"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

failures <- character(0)
check <- function(ok, message) if(!isTRUE(ok)) failures <<- c(failures, message)
payload <- function(widget) jsonlite::fromJSON(as.character(widget$x))

result <- try({
    original_code <- get_genetic_code()
    use_genetic_code("Standard nuclear")
    mk <- function(text, type = "phylogeny") {
        path <- tempfile(fileext = ".newick")
        writeLines(text, path)
        Phylogeny(path, type = type)
    }
    tree <- mk("(A:1.0,B:1.0);")
    modes <- mk("(A:0,B:1);", "mode")

    pi <- F1x4(c(T=.21, C=.29, A=.31, G=.19))
    homogeneous <- GY94(.4, 2.2, pi, scaling_type = "standard")
    branch <- GY94BranchModel(
        list(GY94(.3, 2.2, pi, scaling_type = "standard"),
             GY94(2.5, 2.2, pi, scaling_type = "standard")),
        modes, start_mode = 0)

    labels <- c("slow", "branch")
    assignment <- c("slow", "branch", "slow", "branch", "branch", "slow")
    site_model <- GY94SiteModel(list(homogeneous, branch),
                                site_classes = assignment, labels = labels)
    check(identical(site_model$assignment, match(assignment, labels)),
          "arbitrary character site assignments were reordered")
    check(identical(site_model$n_sites, c(3L, 3L)), "site-class counts are wrong")
    check(inherits(DiscreteGY94(list(homogeneous, branch),
                               site_classes = assignment, labels = labels), "GY94SiteModel"),
          "DiscreteGY94 alias does not construct a GY94SiteModel")

    # Explicit root states avoid making metadata assertions depend on root draws.
    root <- Sequence("TTTGCTAAACCGGACTGG", type = "codon")
    set_palantir_seed(20260825)
    sim <- simulate_gy94_site_model(tree, site_model, sequence = root,
                                    tolerance = 1, rescale_method = "exact")

    expected <- data.frame(
        site = 0:5,
        site_class = assignment,
        class_index = match(assignment, labels) - 1L,
        stringsAsFactors = FALSE)
    check(ncol(sim$alignment) == 6L, "combined alignment length is not six")
    check(isTRUE(all.equal(sim$site_classes, expected, check.attributes = FALSE)),
          "simulation site-class table does not preserve original arbitrary order")
    check(isTRUE(all.equal(attr(sim$alignment, "site_classes"), expected,
                           check.attributes = FALSE)),
          "alignment site-class attribute differs from simulation metadata")
    check(identical(names(sim$site_class_models), labels), "site-class model map has wrong labels")

    if(nrow(sim$substitutions) == 0L) {
        check(FALSE, "inconclusive simulation produced no substitutions")
    } else {
        mapped <- expected$site_class[sim$substitutions$site + 1L]
        mapped_index <- expected$class_index[sim$substitutions$site + 1L]
        check(identical(as.character(sim$substitutions$site_class), mapped),
              "substitution rows carry the wrong site class")
        check(identical(as.integer(sim$substitutions$class_index), mapped_index),
              "substitution rows carry the wrong zero-based class index")
    }

    check(is.list(sim$intervals) && !is.data.frame(sim$intervals) &&
          length(sim$intervals) == 6L && identical(names(sim$intervals), as.character(0:5)),
          "mixed site simulation does not have one named interval entry per site")
    if(is.list(sim$intervals) && length(sim$intervals) == 6L) {
        slow_sites <- which(assignment == "slow")
        branch_sites <- which(assignment == "branch")
        check(all(vapply(sim$intervals[slow_sites], is.null, logical(1))),
              "homogeneous-class sites unexpectedly have branch intervals")
        check(all(vapply(sim$intervals[branch_sites], is.data.frame, logical(1))),
              "branch-class sites lost their interval tables")
    }

    # Alignment subsetting must subset, not drop or renumber, the class metadata.
    selected <- c(5, 2, 6)
    subset <- sim$alignment[, selected]
    subset_meta <- attr(subset, "site_classes")
    check(isTRUE(all.equal(subset_meta, expected[selected, ], check.attributes = FALSE)),
          "Alignment subsetting dropped or misaligned site-class metadata")

    alignment_payload <- payload(plot(sim$alignment))
    check(is.data.frame(alignment_payload$site_classes) &&
          identical(as.character(alignment_payload$site_classes$site_class), assignment),
          "AlignmentPlot payload does not expose site classes")
    phylo_payload <- suppressWarnings(payload(plot(sim, sites = 1)))
    check("site_class" %in% names(phylo_payload$substitutions) &&
          "class_index" %in% names(phylo_payload$substitutions),
          "substitution-history widget payload omits site-class metadata")

    counted <- GY94SiteModel(list(a = homogeneous, b = branch), n_sites = c(2, 3))
    check(identical(counted$assignment, c(1L,1L,2L,2L,2L)) &&
          identical(counted$n_sites, c(2L,3L)), "count-based site assignment is wrong")

    bad_sequence <- Sequence("TTTGCT", type = "codon")
    mismatch <- tryCatch({
        simulate_gy94_site_model(tree, site_model, sequence = bad_sequence)
        FALSE
    }, error = function(e) grepl("one state per simulated site", conditionMessage(e), fixed = TRUE))
    check(mismatch, "wrong-length root sequence was not rejected clearly")

    use_genetic_code(original_code)
    if(length(failures)) fail(paste(failures, collapse = "; "))
    pass()
}, silent = TRUE)

try(use_genetic_code(if(exists("original_code")) original_code else "Standard nuclear"), silent = TRUE)
if(inherits(result, "try-error")) fail(gsub("[\r\n]+", " ", as.character(result)))
