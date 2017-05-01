.globals <- new.env()

.globals$genetic_code_name <- "Standard nuclear"

use_genetic_code <- function(name) {
    .globals$genetic_code_name <- name
}

get_genetic_code <- function() .globals$genetic_code_name

plot.SubstitutionModel <- function(model, prop = 0.25, label_font_size = 12, axis_cex = 1, ...) {

    # Transition rates plot
    tr_plot <- Matrix::image(Matrix::Matrix(model$transition),
                  colorkey = T, useAbs = F,
                  sub = list("Transition rates", fontsize = label_font_size),
                  scales = list(x = list(cex = axis_cex), y = list(cex = axis_cex)),
                  xlab = "", ylab = "", ...)

    # Equilibrium distribution plot
    label_locations <- seq(model$n_states, 1, by = -10) + 1
    label_locations <- label_locations[2:length(label_locations)]

    eq_plot <- lattice::barchart(setNames(rev(model$equilibrium), model$n_states:1),
        scales = list(x = list(cex = axis_cex), y = list(cex = axis_cex, at = label_locations)),
        col = rgb(0.1, 0.4, 0.9), sub = list("Equilibrium", fontsize = label_font_size),
        xlab = "", ylab = "", ...)

    # Combine plots
    print(eq_plot, position = c(0, 0, prop, 1), more = T)
    print(tr_plot, position = c(prop, 0, 1, 1))
}
