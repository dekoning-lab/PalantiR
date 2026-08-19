\name{genetic-code}
\alias{GeneticCode}
\alias{use_genetic_code}
\alias{get_genetic_code}
\title{Genetic code}
\usage{
use_genetic_code(name)

get_genetic_code()
}
\arguments{
\item{name}{Name of a genetic code, for example
\code{"Vertebrate mitochondrial"} or \code{"Standard"}.}
}
\description{
Set and query the genetic code used to build codon models and to translate
codons. The code determines the number of sense codons, which is the dimension
of every codon rate matrix, so it must be set before models are constructed.
}
\details{
Names differ from those used in the 2016 version of this package, which used
\code{"vertebrate_mitochondrial"} in place of \code{"Vertebrate mitochondrial"}.
The vertebrate mitochondrial code has 60 sense codons, not 61.
}
\examples{
use_genetic_code("Vertebrate mitochondrial")
}
