.globals <- new.env()

.globals$genetic_code_name <- "Standard nuclear"

use_genetic_code <- function(name) {
    .globals$genetic_code_name <- name
}

get_genetic_code <- function() .globals$genetic_code_name

Phylogeny <- function(file = NULL, text = NULL) {
    if(!is.null(file)) {
        if(file.exists(file)) {
            text <- scan(file = file, what = character(), quiet = T)
        } else {
            stop("Specified file ", file, " does not exist")
        }
    } else if (is.null(file) && is.null(text)) {
        stop("Either `file` or `text` must be specified")
    }
    PhylogenyClass$new(text)
}
