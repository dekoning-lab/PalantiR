# GY4 -- exact construction of the four Zhang--Nielsen--Yang branch-site classes.

suppressMessages(library(PalantiR))

ID <- "GY4"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

failures <- character(0)
check <- function(ok, message) if(!isTRUE(ok)) failures <<- c(failures, message)
near <- function(x, y, tol = 2e-12) max(abs(as.numeric(x) - as.numeric(y))) < tol

result <- try({
    original_code <- get_genetic_code()
    use_genetic_code("Standard nuclear")
    path <- tempfile(fileext = ".newick")
    writeLines("(background:0,foreground:1);", path)
    modes <- Phylogeny(path, type = "mode")

    background_f <- c(T=.10, C=.20, A=.30, G=.40)
    foreground_f <- rbind(
        c(T=.40, C=.30, A=.20, G=.10),
        c(T=.15, C=.25, A=.35, G=.25),
        c(T=.25, C=.15, A=.20, G=.40))
    counts <- c(2, 3, 4, 5)
    omega0 <- .17
    omega2 <- 4.2
    kappa <- 3.3

    zny <- ZNYBranchSiteModel(
        mode_phylogeny = modes, site_counts = counts,
        omega0 = omega0, omega2 = omega2, kappa = kappa,
        frequencies = background_f, frequency_model = "F1x4",
        foreground_frequencies = foreground_f,
        foreground_frequency_model = "F3x4",
        scaling_type = "synonymous", start_mode = 0)

    check(inherits(zny, "GY94SiteModel"), "ZNY helper did not return a GY94SiteModel")
    check(identical(zny$labels, c("0", "1", "2a", "2b")), "ZNY class labels are wrong")
    check(identical(zny$n_sites, as.integer(counts)), "ZNY site counts are wrong")
    check(identical(zny$assignment, rep(1:4, counts)), "ZNY site assignment is wrong")
    check(all(vapply(zny$models, inherits, logical(1), "GY94BranchModel")),
          "a ZNY site class is not branch heterogeneous")

    expected_omega <- rbind(
        `0`  = c(omega0, omega0),
        `1`  = c(1, 1),
        `2a` = c(omega0, omega2),
        `2b` = c(1, omega2))
    observed_omega <- t(vapply(zny$models, function(spec) {
        vapply(spec$models, `[[`, numeric(1), "omega")
    }, numeric(2)))
    check(near(observed_omega, expected_omega),
          paste("ZNY omega mapping is wrong:", paste(observed_omega, collapse = ",")))

    bg_pi <- F1x4(background_f)
    fg_pi <- F3x4(foreground_f)
    for(label in names(zny$models)) {
        spec <- zny$models[[label]]
        check(spec$start_mode == 0L, paste("wrong start mode in class", label))
        check(near(spec$models[[1]]$equilibrium, bg_pi),
              paste("background frequencies wrong in class", label))
        check(near(spec$models[[2]]$equilibrium, fg_pi),
              paste("foreground frequencies wrong in class", label))
        check(all(vapply(spec$models, `[[`, numeric(1), "kappa") == kappa),
              paste("kappa was not propagated in class", label))
        check(all(vapply(spec$models, `[[`, character(1), "scaling_type") == "synonymous"),
              paste("scaling type was not propagated in class", label))
        check(identical(spec$mode_phylogeny$newick, modes$newick),
              paste("mode phylogeny was not propagated in class", label))
    }

    # The null fixes the foreground positive-selection category at omega = 1.
    null <- ZNYBranchSite(modes, rep(1, 4), omega0 = omega0, omega2 = 1)
    null_omega <- t(vapply(null$models, function(spec) {
        vapply(spec$models, `[[`, numeric(1), "omega")
    }, numeric(2)))
    expected_null <- rbind(c(omega0, omega0), c(1,1), c(omega0,1), c(1,1))
    check(near(null_omega, expected_null), "ZNY branch-site null omega mapping is wrong")

    bad_length <- tryCatch({ ZNYBranchSiteModel(modes, c(1,1,1)); FALSE },
                           error = function(e) grepl("0, 1, 2a, and 2b", conditionMessage(e), fixed = TRUE))
    bad_count <- tryCatch({ ZNYBranchSiteModel(modes, c(1,0,1,1)); FALSE },
                          error = function(e) grepl("positive", conditionMessage(e), fixed = TRUE))
    bad_mode <- tryCatch({
        tree_path <- tempfile(fileext = ".newick")
        writeLines("(A:1,B:1);", tree_path)
        ZNYBranchSiteModel(Phylogeny(tree_path), rep(1,4))
        FALSE
    }, error = function(e) grepl("mode Phylogeny", conditionMessage(e), fixed = TRUE))
    check(bad_length && bad_count && bad_mode, "ZNY helper input validation is incomplete")

    use_genetic_code(original_code)
    if(length(failures)) fail(paste(failures, collapse = "; "))
    pass()
}, silent = TRUE)

try(use_genetic_code(if(exists("original_code")) original_code else "Standard nuclear"), silent = TRUE)
if(inherits(result, "try-error")) fail(gsub("[\r\n]+", " ", as.character(result)))
