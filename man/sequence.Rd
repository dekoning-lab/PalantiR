\name{sequence}
\alias{Sequence}
\alias{sample_sequence}
\alias{as_amino_acid}
\alias{as_compound}
\alias{join}
\title{Sequences and state conversion}
\usage{
sample_sequence(model, length)

as_amino_acid(codons)
}
\arguments{
\item{model}{A \code{SubstitutionModel}; states are drawn from its equilibrium
distribution.}
\item{length}{Number of sites to draw.}
\item{codons}{Codon states, or an \code{Alignment} of them.}
}
\value{
\code{sample_sequence} returns a \code{Sequence}. \code{as_amino_acid} translates
codons to amino acids.
}
\description{
Draw a starting sequence from a model's equilibrium distribution, and convert
between codon and amino-acid representations.
}
\details{
\code{as_amino_acid} operates elementwise and returns a plain character vector,
discarding the matrix dimensions and taxon names of an \code{Alignment}. To
translate an alignment while keeping its shape, assign back into a copy:
\code{aa <- aln; aa[] <- as_amino_acid(aln)}.
}
\seealso{\code{\link{as.fasta}}}
