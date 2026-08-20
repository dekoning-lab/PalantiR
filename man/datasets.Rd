\name{datasets}
\alias{amino_acid_equilibrium}
\alias{amino_acid_fitness_N_1000}
\alias{psi_c50_1}
\title{Example amino-acid distributions and fitnesses}
\format{
\code{amino_acid_equilibrium}: a 6 x 20 numeric matrix of amino-acid stationary
distributions, one site per row, columns named for the amino acids in this
package's order \code{ACDEFGHIKLMNPQRSTVWY}. Rows sum to 1.

\code{amino_acid_fitness_N_1000}: a 6 x 20 numeric matrix of relative fitnesses
for the same six sites, same column order, each row scaled so that its most
frequent amino acid has fitness 1.

\code{psi_c50_1}: a 50 x 20 numeric matrix of amino-acid distributions, one per
row, with no dimnames. Rows sum to 1, so these are frequency profiles rather
than fitnesses; pass them through \code{\link{equilibrium_to_fitness}} before
building a \code{\link{MutationSelection}} model.
}
\description{
Example inputs shipped with the package. \code{amino_acid_equilibrium} holds
amino-acid stationary distributions, \code{amino_acid_fitness_N_1000} the
corresponding relative fitnesses at an effective population size of 1000, and
\code{psi_c50_1} a set of 50 amino-acid profiles.
}
\details{
The raw files these were built from are in \code{inst/extdata}
(\code{amino_acid_equilibrium.csv}, \code{amino_acid_fitness_N_1000.csv} and
\code{psi-c50-1.txt}), reached with \code{system.file("extdata",
"amino_acid_equilibrium.csv", package = "PalantiR")}. Every line of
\code{psi-c50-1.txt} ends in a tab, which \code{read.delim} turns into a 21st,
all-zero column; the packaged \code{psi_c50_1} has that column removed.

Example trees are in \code{inst/extdata} too, and are reached with
\code{system.file("extdata", "primates.newick", package = "PalantiR")}.
}
\note{
Before version 1.2.0 these datasets shipped as \code{.csv} and \code{.txt} files
in \code{data/}, where \code{data()} read them with \code{read.table}'s
\code{sep = ";"} convention and delivered unusable single-column objects. They
are now \code{.rda} files, and the dataset once documented as \code{psi-c50-1}
has been renamed \code{psi_c50_1} so that it can be used as an R name.
}
\examples{
data(amino_acid_fitness_N_1000)
dim(amino_acid_fitness_N_1000)
}
\keyword{datasets}
