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
                                 segment_length = 0.001, tolerance = 0.001,
                                 rescale_method = "segments")

simulate_with_shared_substitution_heterogeneity(phylogeny, switching_model,
                                 substitution_models, sequence, start_mode,
                                 rate = 1, switching_rate = 1,
                                 segment_length = 0.0001, tolerance = 0.001)

simulate_with_shared_time_heterogeneity(phylogeny, switching_model,
                                 substitution_models, sequence, start_mode,
                                 rate = 1, switching_rate = 1,
                                 segment_length = 0.0001, tolerance = 0.001)
}
\arguments{
\item{phylogeny}{A \code{Phylogeny} of type \code{"phylogeny"} carrying the
branch lengths to simulate over.}
\item{model}{A single \code{SubstitutionModel}.}
\item{models, substitution_models}{List of models, one per mode. See Details for
indexing. All of them must share a \code{type} and a \code{scaling_type}.}
\item{mode_phylogeny}{A \code{Phylogeny} of the same topology whose branch
"lengths" are integer mode labels, assigning each branch to a model. It must be
built with \code{Phylogeny(path, type = "mode")}; a branch-length tree is now
rejected rather than silently run under \code{models[[1]]}.}
\item{switching_model}{A \code{SubstitutionModel} over the modes themselves,
whose states index \code{substitution_models}. The mode process runs along the
tree and drives the changes of model.}
\item{sequence}{Starting sequence, from \code{\link{sample_sequence}}. Its length
is the number of independent sites simulated.}
\item{start_mode}{Mode in force at the root.}
\item{rate}{Scalar rate multiplier applied to every branch length. Use this for
among-site rate variation, one call per site.}
\item{switching_rate}{Rate multiplier for the mode process, applied to the
branch lengths seen by \code{switching_model} only.}
\item{segment_length}{Step size for the transient renormalisation.}
\item{tolerance}{Convergence criterion for the transient renormalisation.}
\item{rescale_method}{Which renormalisation to use on a branch whose mode has
just changed, either \code{"segments"} (the default) or \code{"exact"}. See
Details.}
}
\value{
A list of class \code{Simulation}. \code{simulate_over_phylogeny} returns
elements \code{phylogeny}, \code{model}, \code{substitutions}, \code{alignment},
\code{intervals} and \code{type}; the three heterogeneous simulators return the
list of models under \code{models} rather than a single \code{model}.

\code{intervals} is \code{NULL} for \code{simulate_over_phylogeny}. For
\code{simulate_over_interval_phylogeny} and
\code{simulate_with_shared_substitution_heterogeneity} it is a single data frame
of mode intervals, shared by every site. For
\code{simulate_with_shared_time_heterogeneity}, where the mode process runs
independently at each site, it is a \emph{list} of such data frames, one per
site.

The \code{substitutions} data frame records every substitution: the site, the
node below which it occurred, its time on the branch, the codons and amino acids
involved, and whether it was synonymous. Its \code{site} column is
\strong{0-based}: the first site of the sequence is site \code{0}. This record
cannot be recovered from the alignment afterwards, so retain it if any
downstream analysis needs substitution histories.
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

\strong{Where the modes come from.} \code{simulate_over_interval_phylogeny} takes
them from a guide tree, so they are fixed in advance and shared by every site.
The two \code{simulate_with_shared_*} simulators instead sample them from
\code{switching_model}: \code{shared_substitution_heterogeneity} samples one mode
history for the whole tree, which all sites then share, while
\code{shared_time_heterogeneity} samples an independent mode history for each
site.

\strong{Transient renormalisation.} On a branch whose mode differs from its
parent's, the sequence begins out of equilibrium with respect to its new model,
and the realised substitution rate departs from the model's equilibrium rate. Two
renormalisations are available, chosen with \code{rescale_method}.

\code{"segments"}, the default, splits the branch into steps of
\code{segment_length}, advances the transient distribution along each, and
renormalises the rate matrix against that distribution rather than the
equilibrium one, so that branch lengths retain their meaning across the change.
It disengages once the transient distribution reaches the mode's equilibrium
within \code{tolerance}. Setting \code{tolerance} to 1 or more keeps it out of
the path entirely, giving an instantaneous switch. Its accuracy and its cost are
both set by \code{segment_length}.

\code{"exact"} does the same job without discretising. The branch's budget of
scaled-class events is the one the equilibrium model would deliver over the
branch length; the transient model is run for however much intrinsic time it
takes to deliver exactly that budget, solving the time change
\eqn{\tau \mapsto \int_0^{\tau} f(u)\,du} for the intrinsic duration, where
\eqn{f} is the scaled-class flux under the transient distribution. Event times
are then mapped back onto the branch. This is exact in the segment size and
usually much faster, but it is defined only when the models carry a scaled class:
with \code{scaling_type = "none"} it now raises an error instead of quietly
running the segment scheme. It does not use \code{segment_length}, which
therefore has no effect on a run made with \code{rescale_method = "exact"}.

Either renormalisation is defined for single-codon models only. A
\code{CoEvolution} model reaching it raises an error naming the mismatch; pass a
\code{tolerance} of 1 or more to keep it out of the path.

Only \code{simulate_over_interval_phylogeny} exposes \code{rescale_method}; the
\code{simulate_with_shared_*} simulators always use \code{"segments"}, and their
much smaller default \code{segment_length} makes them correspondingly expensive
on long branches.

Note that \code{tolerance} is an absolute root-mean-square deviation on a
distribution whose entries scale as one over the number of states, so its
effective strictness depends on the size of the state space.

\strong{Reproducibility.} The simulation engine has its own random number
generator, which \code{set.seed} does not reach. Use
\code{\link{set_palantir_seed}}.
}
\seealso{\code{\link{MutationSelection}}, \code{\link{Phylogeny}},
\code{\link{set_palantir_seed}}}
