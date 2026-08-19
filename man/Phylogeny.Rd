\name{Phylogeny}
\alias{Phylogeny}
\alias{phylogeny_to_intervals}
\alias{compare_modes}
\alias{plot.Phylogeny}
\alias{print.Phylogeny}
\title{Read a phylogeny}
\usage{
Phylogeny(newick_path, type = "phylogeny")
}
\arguments{
\item{newick_path}{Path to a file containing a Newick string. Note this is a
\emph{file path}, not the string itself.}
\item{type}{Either \code{"phylogeny"} or \code{"mode"}.}
}
\value{
A list of class \code{Phylogeny} with elements \code{json} (the tree as nested
JSON, carrying node indices, labels and branch lengths), \code{string},
\code{newick}, \code{n_nodes} and \code{type}. Node indices in \code{json} are
the same indices reported in a simulation's \code{substitutions} table, so the
JSON is the way to map substitutions back to named branches.
}
\description{
Read a phylogeny from a Newick file. A guide tree assigning branches to
substitution models is read the same way, with integer mode labels in place of
branch lengths.
}
\examples{
\dontrun{
ph <- Phylogeny("tree.newick")
ph
str(jsonlite::fromJSON(ph$json), max.level = 2)
}
}
