\name{equilibrium_to_fitness}
\alias{equilibrium_to_fitness}
\title{Convert a stationary distribution to relative fitnesses}
\usage{
equilibrium_to_fitness(equilibrium, population_size,
                       nucleotide_equilibrium = NULL)
}
\arguments{
\item{equilibrium}{Amino-acid stationary distribution, length 20, in this
package's order \code{ACDEFGHIKLMNPQRSTVWY}.}
\item{population_size}{The effective population size at which to interpret the
distribution.}
\item{nucleotide_equilibrium}{Optionally the nucleotide equilibrium of the
mutation process, length 4, in the order \code{T}, \code{C}, \code{A},
\code{G} and strictly positive. Supplying it makes the inversion
degeneracy-aware; see Details.}
}
\value{
A 20 x 1 numeric matrix of relative fitnesses, the fittest amino acid taking
fitness exactly 1. It is a one-column matrix rather than a plain vector, which
\code{\link{MutationSelection}} accepts as its \code{fitness} argument;
\code{as.numeric} it if a vector is wanted.
}
\description{
Recover relative amino-acid fitnesses from a stationary distribution under the
mutation-selection model of Halpern and Bruno (1998), by inverting
\eqn{\pi_b p_{ba} / (\pi_a p_{ab}) = e^{2 N_e s}}.
}
\details{
This is a \emph{de-scaling} step. What an alignment constrains is the
population-scaled fitness \eqn{2 N_e s}, so a population size must be supplied to
recover \eqn{s} itself. If the intention is to hold the fitness landscape fixed
while varying \eqn{N_e} across a tree, supply the branch-length-weighted mean
\eqn{N_e}; that conserves the tree-averaged population-scaled fitness between a
time-homogeneous simulation and its time-heterogeneous counterpart.

The most frequent amino acid is taken as the reference and assigned fitness 1.

\strong{The degeneracy tilt, and \code{nucleotide_equilibrium}.} With
\code{nucleotide_equilibrium} left at \code{NULL} the inversion is the plain one,
\eqn{\psi_A \propto \log \pi_A / 2 N_e}, which ignores the mutational weight of
each amino acid's codons. Amino acids with more (or more
mutationally-favoured) codons then come out over-represented: a
\code{MutationSelection} model built on the result has an amino-acid marginal
that is the degeneracy-tilted transform of \code{equilibrium}, not
\code{equilibrium} itself. Under a uniform target and a uniform nucleotide
equilibrium the realised marginal spreads by exactly a factor of six, the ratio
of the six-codon to the one-codon families.

Supplying \code{nucleotide_equilibrium} divides that weight out, so that a
\code{MutationSelection} model built on the returned fitnesses \emph{with the
same nucleotide equilibrium} reproduces the requested amino-acid marginal. The
default is unchanged for backward compatibility: an existing pipeline is
internally consistent under the plain convention, and switching it changes the
fitted landscape.
}
\references{
Halpern, A. L. and Bruno, W. J. (1998) Evolutionary distances for protein-coding
sequences: modeling site-specific residue frequencies.
\emph{Molecular Biology and Evolution} \bold{15}, 910--917.
}
\seealso{\code{\link{MutationSelection}}}
