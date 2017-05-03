.globals <- new.env()

.globals$genetic_code_name <- "Standard nuclear"

use_genetic_code <- function(name) {
    .globals$genetic_code_name <- name
}

get_genetic_code <- function() .globals$genetic_code_name
