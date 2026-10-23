\name{datasets}
\alias{amino_acid_equilibrium}
\alias{amino_acid_fitness_N_1000}
\alias{psi_c50_1}
\alias{psi_c50_2}
\title{Example amino-acid distributions and fitnesses}
\format{
\code{amino_acid_equilibrium}: a 6 x 20 numeric matrix of amino-acid stationary
distributions, one site per row, columns named for the amino acids in this
package's order \code{ACDEFGHIKLMNPQRSTVWY}. Rows sum to 1.

\code{amino_acid_fitness_N_1000}: a 6 x 20 numeric matrix of relative fitnesses
for the same six sites, same column order, each row scaled so that its most
frequent amino acid has fitness 1.

\code{psi_c50_1}, \code{psi_c50_2}: 50 x 20 numeric matrices of amino-acid
frequency profiles, one per row, columns named in the same order and no row
names; row \eqn{k} is profile \eqn{k} of the source table. Rows sum to 1, so
these are frequency profiles rather than fitnesses; pass them through
\code{\link{equilibrium_to_fitness}} before building a
\code{\link{MutationSelection}} model.
}
\description{
Example inputs shipped with the package. \code{amino_acid_equilibrium} holds
amino-acid stationary distributions, \code{amino_acid_fitness_N_1000} the
corresponding relative fitnesses at an effective population size of 1000, and
\code{psi_c50_1} and \code{psi_c50_2} two independent estimates of the 50
amino-acid profiles of the C50 profile mixture of Le, Gascuel and Lartillot
(2008).
}
\details{
\code{psi_c50_1} is Table 4 (\dQuote{first run}) and \code{psi_c50_2} is
Table 5 (\dQuote{second run}) of the supplementary material of Le, Gascuel and
Lartillot (2008). The tables print each profile as whole percentages, so a
printed row sums to between 95 and 103; each row here is the printed row
divided by its own total, and multiplying it by that total recovers the printed
integers. Table 5 is ordered to match Table 4, and the paper identifies its
first 41 profiles unequivocally with those of the first run. The tables give no
mixture weights, so neither data set carries any.

The raw files these were built from are in \code{inst/extdata}
(\code{amino_acid_equilibrium.csv}, \code{amino_acid_fitness_N_1000.csv},
\code{psi-c50-1.txt} and \code{psi-c50-2.txt}), reached with
\code{system.file("extdata", "amino_acid_equilibrium.csv", package =
"PalantiR")}. The two \code{psi-c50} files are tab-separated, with no header
and one profile per line. Each line has 21 fields: the 20 frequencies in
\code{ACDEFGHIKLMNPQRSTVWY} order, then a final \code{0}, the stop-codon
symbol \code{Z} of Palantir-Perl, which read these files as
\code{ACDEFGHIKLMNPQRSTVWYZ}. The packaged matrices drop that field.

Example trees are in \code{inst/extdata} too, and are reached with
\code{system.file("extdata", "primates.newick", package = "PalantiR")}.
}
\source{
Le, S. Q., Gascuel, O. and Lartillot, N. (2008) Empirical profile mixture
models for phylogenetic reconstruction. \emph{Bioinformatics} \bold{24},
2317--2323. \doi{10.1093/bioinformatics/btn445}. Supplementary material
(\file{bioinf-2008-0768-File009.pdf}), Table 4 (\code{psi_c50_1}) and Table 5
(\code{psi_c50_2}).
}
\note{
Before version 1.3.0, \code{psi_c50_1} had its V and W columns exchanged:
column 18, documented as V, held Table 4's W values, and column 19, documented
as W, held its V values. Profile 2, which the paper gives as 41\% I and 43\% V,
therefore read as 41\% I and 43\% W. Every other column matched the table. The
swap dates to the file's first appearance, as \file{examples/psi-c50-1.txt} in
Palantir-Perl (2015), from which PalantiR copied it unchanged; Palantir-Perl's
\file{examples/psi-c50-2.txt}, the source of \code{psi_c50_2}, carries the same
swap. Both raw files and both matrices now follow the tables, and every row of
each reproduces its table exactly. Simulations run with \code{psi_c50_1} under
earlier versions used profiles with V and W exchanged.

That the tables, and not the files, label V and W correctly is supported by
composition rather than by independent confirmation. With the tables' labels,
V averages 6.8\% and W 0.6\% across the profiles of Table 4, typical of
proteins, and across profiles V rises and falls with I (correlation 0.88) and W
with the other aromatics, F and Y (0.74 each), as chemically similar residues
should. With the files'
labels these relationships vanish, and W reaches 43\% in an isoleucine profile.

Before version 1.2.0 these datasets shipped as \code{.csv} and \code{.txt} files
in \code{data/}, where \code{data()} read them with \code{read.table}'s
\code{sep = ";"} convention and delivered unusable single-column objects. They
are now \code{.rda} files, and the dataset once documented as \code{psi-c50-1}
has been renamed \code{psi_c50_1} so that it can be used as an R name.
}
\examples{
data(amino_acid_fitness_N_1000)
dim(amino_acid_fitness_N_1000)

# Profile 2 of Table 4: 41\% I, 43\% V, no W
data(psi_c50_1)
round(100 * psi_c50_1[2, c("I", "V", "W")], 1)
}
\keyword{datasets}
