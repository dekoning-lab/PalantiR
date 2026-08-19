\name{equilibrium_to_fitness}
\alias{equilibrium_to_fitness}
\title{Convert a stationary distribution to relative fitnesses}
\usage{
equilibrium_to_fitness(equilibrium, population_size)
}
\arguments{
\item{equilibrium}{Amino-acid stationary distribution, length 20, in this
package's order \code{ACDEFGHIKLMNPQRSTVWY}.}
\item{population_size}{The effective population size at which to interpret the
distribution.}
}
\value{Numeric vector of 20 relative fitnesses.}
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
}
\references{
Halpern, A. L. and Bruno, W. J. (1998) Evolutionary distances for protein-coding
sequences: modeling site-specific residue frequencies.
\emph{Molecular Biology and Evolution} \bold{15}, 910--917.
}
\seealso{\code{\link{MutationSelection}}}
