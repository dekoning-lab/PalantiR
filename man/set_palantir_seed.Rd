\name{set_palantir_seed}
\alias{set_palantir_seed}
\title{Seed the simulation engine}
\usage{
set_palantir_seed(seed)
}
\arguments{
\item{seed}{A whole number between 0 and 4294967295.}
}
\value{
Invisibly \code{NULL}, called for its effect on the engine.
}
\description{
Make a run of PalantiR simulations reproducible.
}
\details{
The simulation engine draws from its own C++ Mersenne twister, not from R's
random number generator. \strong{\code{set.seed} therefore has no effect on any
simulation}: waiting times, jumps, sampled sequences and mode histories are all
untouched by it. Call \code{set_palantir_seed} instead, once, before the
simulations to be reproduced.

A session that never calls it seeds the engine from \code{std::random_device},
so unseeded runs stay non-deterministic. Everything drawn after a given call to
\code{set_palantir_seed} is a deterministic function of that seed and of the
order in which the draws are made, so reproducing a result means reproducing the
whole sequence of calls that led to it, not only the last one.

Before version 1.2.0 there was no way to seed the engine at all, and each
translation unit held a separately seeded copy of the generator; no PalantiR
simulation was reproducible.
}
\examples{
set_palantir_seed(42)
}
\seealso{\code{\link{simulate_over_phylogeny}}}
