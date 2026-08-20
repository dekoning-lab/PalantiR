\name{genetic-code}
\alias{GeneticCode}
\alias{use_genetic_code}
\alias{get_genetic_code}
\alias{genetic_code_names}
\title{Genetic code}
\usage{
use_genetic_code(name)

get_genetic_code()

genetic_code_names()
}
\arguments{
\item{name}{Name of a genetic code, for example
\code{"Vertebrate mitochondrial"} or \code{"Standard nuclear"}. It must be one of
the names returned by \code{genetic_code_names()}.}
}
\value{
\code{genetic_code_names} returns a character vector of the valid code names,
read from the package's own code table. \code{get_genetic_code} returns the name
of the code currently in force, \code{"Standard nuclear"} by default.
}
\description{
Set and query the genetic code used to build codon models and to translate
codons. The code determines the number of sense codons, which is the dimension
of every codon rate matrix, so it must be set before models are constructed.
}
\details{
\code{use_genetic_code} checks its argument against \code{genetic_code_names()}
and errors on anything else. It used to accept any value at all, including
numbers and \code{NULL}, and the mistake surfaced much later as an unrelated
error from the first model built afterwards.

Nothing binds a model to the code it was built under. Changing the code after
building models leaves those models intact but reinterprets their state indices,
which mistranslates codons or fails with "Unexpected integer in codon sequence".
Set the code once, at the top of the script.

Names differ from those used in the 2016 version of this package, which used
\code{"vertebrate_mitochondrial"} in place of \code{"Vertebrate mitochondrial"}.
The vertebrate mitochondrial code has 60 sense codons, not 61.
}
\examples{
genetic_code_names()
use_genetic_code("Vertebrate mitochondrial")
get_genetic_code()
use_genetic_code("Standard nuclear")
}
