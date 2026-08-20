# Runs every regression test in tests/review/. Each test file is a stand-alone
# script that prints "<ID>: PASS ..." or "<ID>: FAIL ..." and always exits 0,
# so each is run in a subprocess and judged by its output. The suite fails
# (and with it R CMD check) if any test reports FAIL, reports nothing, or
# exits abnormally.

files <- sort(list.files(file.path("review"), pattern = "\\.R$", full.names = TRUE))
stopifnot(length(files) > 0)

rscript <- file.path(R.home("bin"), "Rscript")
libs <- paste(.libPaths(), collapse = .Platform$path.sep)

failures <- character(0)
for (f in files) {
    id <- sub("\\.R$", "", basename(f))
    t0 <- proc.time()[3]
    out <- suppressWarnings(system2(rscript, shQuote(f), stdout = TRUE, stderr = TRUE,
                                    env = paste0("R_LIBS=", libs)))
    status <- attr(out, "status"); if (is.null(status)) status <- 0
    text <- paste(out, collapse = "\n")
    skipped <- grepl(paste0(id, ": SKIP"), text)
    ok <- status == 0 && !grepl("FAIL", text) &&
          (grepl(paste0(id, ": PASS"), text) || skipped)
    cat(sprintf("%-6s %s  (%.1fs)\n", id,
        if (!ok) "FAIL" else if (skipped) "SKIP" else "PASS", proc.time()[3] - t0))
    if (!ok) {
        failures <- c(failures, id)
        cat(text, "\n")
    }
}
cat(sprintf("review suite: %d/%d passed\n", length(files) - length(failures), length(files)))
if (length(failures)) stop("failing tests: ", paste(failures, collapse = ", "))
