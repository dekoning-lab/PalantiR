\name{simulate}
\alias{simulate_over_phylogeny}
\alias{simulate_over_interval_phylogeny}
\alias{simulate_with_shared_substitution_heterogeneity}
\alias{simulate_with_shared_time_heterogeneity}
\alias{print.Simulation}
\title{Simulate sequence evolution along a phylogeny}
\usage{
simulate_over_phylogeny(phylogeny, model, sequence, rate = 1)

simulate_over_interval_phylogeny(phylogeny, mode_phylogeny, models, sequence,
                                 start_mode, rate = 1,
                                 segment_length = 0.001, tolerance = 0.001)
}
\arguments{
\item{phylogeny}{A \code{Phylogeny} carrying the branch lengths to simulate over.}
\item{model}{A single \code{SubstitutionModel}.}
\item{models}{List of models, one per mode. See Details for indexing.}
\item{mode_phylogeny}{A \code{Phylogeny} of the same topology whose branch
"lengths" are integer mode labels, assigning each branch to a model.}
\item{sequence}{Starting sequence, from \code{\link{sample_sequence}}. Its length
is the number of independent sites simulated.}
\item{start_mode}{Mode in force at the root.}
\item{rate}{Scalar rate multiplier applied to every branch length. Use this for
among-site rate variation, one call per site.}
\item{segment_length}{Step size for the transient renormalisation.}
\item{tolerance}{Convergence criterion for the transient renormalisation.}
}
\value{
A list of class \code{Simulation} with elements \code{phylogeny}, \code{model},
\code{substitutions}, \code{alignment}, \code{intervals} and \code{type}. The
\code{substitutions} data frame records every substitution: the node below which
it occurred, its time on the branch, the codons and amino acids involved, and
whether it was synonymous. This record cannot be recovered from the alignment
afterwards, so retain it if any downstream analysis needs substitution histories.
}
\description{
Simulate codon sequence evolution along a phylogeny, either under a single model
or under models that change across the tree.
}
\details{
\strong{Mode indexing.} A mode labelled \eqn{m} on the guide tree selects
\code{models[[m + 1]]}. A guide tree using labels 1 and 2 therefore needs a
three-element list, with an unused first element. Getting this wrong selects the
wrong model without warning.

\strong{Transient renormalisation.} On a branch whose mode differs from its
parent's, the sequence begins out of equilibrium with respect to its new model,
and the realised substitution rate departs from the model's equilibrium rate. The
branch is therefore split into segments of \code{segment_length}, the transient
distribution advanced along each, and the rate matrix renormalised against that
distribution rather than the equilibrium one, so that branch lengths retain their
meaning across the change. It disengages once the transient distribution reaches
the mode's equilibrium within \code{tolerance}. Setting \code{tolerance} to 1 or
more keeps it out of the path entirely, giving an instantaneous switch.

Note that \code{tolerance} is an absolute root-mean-square deviation on a
distribution whose entries scale as one over the number of states, so its
effective strictness depends on the size of the state space.
}
\seealso{\code{\link{MutationSelection}}, \code{\link{Phylogeny}}}
