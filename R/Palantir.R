.globals <- new.env()

.globals$genetic_code_name <- "Standard nuclear"

use_genetic_code <- function(name) {
    .globals$genetic_code_name <- name
}

get_genetic_code <- function() .globals$genetic_code_name

.all_same <- function(x) all(x == x[[1]])
.get_member <- function(x, a) lapply(x, function(x) x[[a]])
.get_attr <- function(x, a) lapply(x, function(x) attr(x, a))

.pad_with_NA <- function(df, col_names) {
    for(name in col_names) {
        if(is.na(match(name, names(df)))) {
            df[[name]] <- NA
        }
    }
    df
}

# Join simulations
join <- function(...) {
    sims <- list(...)
    first <- sims[[1]]
    # check compatibiilty
    alignments <- .get_member(sims, "alignment")
    phylogenies <- .get_member(sims, "phylogeny")
    subs <- .get_member(sims, "substitutions")

    stopifnot(.all_same(.get_member(phylogenies, "newick")))
    stopifnot(.all_same(.get_attr(alignments, "type")))
    stopifnot(.all_same(.get_member(sims, "type")))

    joined_alignment <- do.call(cbind, alignments)
    alignment_type <- attr(first$alignment, "type")
    attr(joined_alignment, "type") <- alignment_type
    class(joined_alignment) <- "Alignment"

    col_names <- Reduce(union, .get_attr(subs, "names"))
    padded <- lapply(subs, .pad_with_NA, col_names)
    joined_subs <- do.call(rbind, padded)

    joined_sim <- list(
        phylogeny = first$phylogeny,
        models = .get_member(sims, "model"),
        substitutions = joined_subs,
        alignment = joined_alignment,
        intervals = .get_member(sims, "intervals"),
        type = first$type
    )
    class(joined_sim) <- "Simulation"
    return(joined_sim)
}

