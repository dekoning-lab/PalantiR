# M7 (R half) -- use_genetic_code() validates its argument at the setter.
#
# Before the fix use_genetic_code() stored whatever it was given. A misspelt
# name, a number or NULL was accepted silently and then poisoned every later
# call: the failure surfaced from inside a model constructor as "Unknown genetic
# code", or from a simulation as "Unexpected integer in codon sequence", with no
# hint that the setter was at fault and no list of the valid names anywhere.
#
# EXPECTED against the OLD library: FAIL -- use_genetic_code("Standard") returns
#   quietly and leaves the package holding an unusable code name.
# EXPECTED against the FIXED library: PASS -- invalid names, non-character
#   values, NULL and length != 1 all error immediately; the message lists the
#   valid names; the active code is left untouched by a rejected call; and every
#   name from genetic_code_names() is still accepted.

.libPaths(c('/Users/jasondk/bot-workspace/incidental-convergence/rlib_review', .libPaths()))
suppressMessages(library(PalantiR))

ID <- "M7R"
done <- function(msg) { cat(msg, "\n", sep = ""); quit(save = "no", status = 0) }
fail <- function(reason) done(sprintf("%s: FAIL %s", ID, reason))
pass <- function() done(sprintf("%s: PASS", ID))

# returns the error message, or NULL if the call succeeded
rejected <- function(expr) {
    tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
}

original <- get_genetic_code()

result <- try({
    valid <- genetic_code_names()

    # 1. a plausible-but-wrong name must be rejected, naming the alternatives
    msg <- rejected(use_genetic_code("Standard"))
    if (is.null(msg)) {
        fail("use_genetic_code(\"Standard\") was accepted")
    }
    if (!grepl("Standard nuclear", msg, fixed = TRUE)) {
        fail(sprintf("the error does not list the valid names: %s", msg))
    }
    if (!grepl("Standard", msg, fixed = TRUE)) {
        fail(sprintf("the error does not name the offending value: %s", msg))
    }
    if (!identical(get_genetic_code(), original)) {
        fail(sprintf("a rejected call still changed the active code to %s",
                     get_genetic_code()))
    }

    # 2. the other ways of getting it wrong
    for (bad in list(42, NULL, NA_character_, c("Standard nuclear", "Standard nuclear"),
                     character(0), list("Standard nuclear"))) {
        if (is.null(rejected(use_genetic_code(bad)))) {
            fail(sprintf("use_genetic_code(%s) was accepted",
                         paste(utils::capture.output(str(bad)), collapse = " ")))
        }
        if (!identical(get_genetic_code(), original)) {
            fail("a rejected call still changed the active code")
        }
    }

    # 3. every advertised name still works, and still builds a code
    for (name in valid) {
        if (!is.null(rejected(use_genetic_code(name)))) {
            fail(sprintf("the valid name \"%s\" was rejected", name))
        }
        if (!identical(get_genetic_code(), name)) {
            fail(sprintf("use_genetic_code(\"%s\") did not take effect", name))
        }
        if (length(GeneticCode()) < 60) {
            fail(sprintf("\"%s\" did not build a usable genetic code", name))
        }
    }

    use_genetic_code(original)
    pass()
}, silent = TRUE)

try(use_genetic_code(original), silent = TRUE)

if (inherits(result, "try-error")) {
    fail(gsub("[\r\n]+", " ", as.character(result)))
}
