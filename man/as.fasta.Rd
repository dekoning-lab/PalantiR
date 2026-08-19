\name{as.fasta}
\alias{as.fasta}
\alias{as.fasta.Alignment}
\title{Write an alignment in FASTA format}
\usage{
as.fasta(x, file = "", ...)
}
\arguments{
\item{x}{An \code{Alignment}, a taxa-by-sites character matrix whose row names
are the taxon names.}
\item{file}{Destination path. The default writes to the console.}
\item{...}{Unused, for consistency with the generic.}
}
\description{
Write an alignment as FASTA, one record per taxon.
}
\examples{
\dontrun{
as.fasta(sim$alignment, file = "codons.fasta")
aa <- sim$alignment
aa[] <- as_amino_acid(sim$alignment)
as.fasta(aa, file = "amino_acids.fasta")
}
}
