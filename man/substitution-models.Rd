\name{substitution-models}
\alias{MutationSelection}
\alias{CoEvolution}
\alias{HasegawaKishinoYano}
\alias{GeneralTimeReversible}
\alias{MarkovModulatedMutationSelection}
\alias{print.SubstitutionModel}
\title{Substitution models}
\usage{
HasegawaKishinoYano(equilibrium, transition_rate = 1, transversion_rate = 1)

GeneralTimeReversible(equilibrium, exchangeability)

MutationSelection(population_size, mutation_rate, nucleotide_model, fitness,
                  scaling_type = "synonymous")

CoEvolution(population_size, mutation_rate, nucleotide_model, fitness_1,
            fitness_2, delta, scaling_type = "synonymous")

MarkovModulatedMutationSelection(mutation_selection_models, switching_model)
}
\arguments{
\item{equilibrium}{Numeric vector of equilibrium frequencies. Length 4 for a
nucleotide model.}
\item{transition_rate, transversion_rate}{Relative rates of transitions and
transversions. Only their ratio matters when the model is used as the mutation
process of a \code{MutationSelection} or \code{CoEvolution} model, since the
codon model is renormalised afterwards. Used directly to simulate nucleotide
evolution the model is not renormalised, so their absolute size sets the overall
rate and scales the meaning of the branch lengths with it.}
\item{exchangeability}{Symmetric matrix of relative exchange rates.}
\item{population_size}{Effective population size, \eqn{N_e}.}
\item{mutation_rate}{Scalar multiplier on the mutation process.}
\item{nucleotide_model}{A nucleotide substitution model, typically from
\code{HasegawaKishinoYano}, supplying the mutation process.}
\item{fitness, fitness_1, fitness_2}{Numeric vector of 20 relative amino-acid
fitnesses. \strong{Order matters and is not checked}: see Details.}
\item{delta}{Matrix of pairwise interaction terms for the co-evolution model.}
\item{scaling_type}{How the rate matrix is normalised. One of
\code{"synonymous"} (the dS gauge), \code{"synonymous-per-codon"},
\code{"non-synonymous"}, \code{"substitution"}, or \code{"none"}.
\code{"dS"} and \code{"ds"} are aliases of \code{"synonymous"};
\code{"standard"} is an alias of \code{"substitution"}. See Details.}
\item{switching_model}{An exchangeable model governing switching between modes.}
\item{mutation_selection_models}{List of mutation-selection models, one per mode.}
}
\value{
A list of class \code{SubstitutionModel} with elements \code{equilibrium},
\code{transition} (the normalised rate matrix), \code{sampling} (its cumulative
form, used by the simulators), \code{scaling} (the normalising constant applied),
\code{scaling_type}, \code{n_states} and \code{type}. dS-scaled models also
report \code{neutral_total_rate}, \code{neutral_synonymous_rate}, and the
Q-weighted \code{synonymous_opportunities} used to define the gauge; these
fields are \code{NA} for the other scaling types. They also report the
per-state neutral-reference \code{dS_outflux} used by the transient rescaler.
}
\description{
Construct substitution models. \code{MutationSelection} implements the
mutation-selection formulation of Halpern and Bruno (1998) over codons, in which
amino-acid fitnesses and an effective population size determine fixation
probabilities and hence substitution rates.
}
\details{
\strong{Amino-acid ordering.} \code{fitness} is interpreted in this package's
own order, \code{ACDEFGHIKLMNPQRSTVWY} (alphabetical by one-letter code). Other
software uses different conventions; PAML and its rate files, for instance, use
\code{ARNDCQEGHILKMFPSTWYV}. The vector is unnamed and its ordering cannot be
validated, so a mismatch produces a silently permuted fitness landscape rather
than an error. Reorder external input explicitly.

\strong{Scaling.} \code{"synonymous"} is the default dS gauge and is stored
canonically as \code{"dS"}. The selected generator is divided by one third of
the total rate of its mutation-only neutral reference. One unit of branch
length is therefore one neutral synonymous substitution per synonymous
opportunity, equivalently one neutral nucleotide substitution per site. The
reference sets amino-acid fitness differences to zero and depends only on the
nucleotide mutation model, so site-heterogeneous fitness profiles sharing that
model also share one clock.

\code{"synonymous-per-codon"} retains the earlier convention in which one
branch-length unit gives one expected synonymous event per codon at the
model's equilibrium. \code{"non-synonymous"} and \code{"substitution"}
similarly count amino-acid-changing events or all events per codon;
\code{"none"} applies no scaling. Choose the convention matching the supplied
tree. An unrecognised name raises an error. \code{"dS"} and \code{"ds"} are
explicit aliases for the dS gauge. For compatibility with the 2016 version of
this package, \code{"standard"} is accepted as an alias and stored canonically
as \code{"substitution"}.
}
\references{
Halpern, A. L. and Bruno, W. J. (1998) Evolutionary distances for protein-coding
sequences: modeling site-specific residue frequencies.
\emph{Molecular Biology and Evolution} \bold{15}, 910--917.

Hasegawa, M., Kishino, H. and Yano, T. (1985) Dating of the human-ape splitting
by a molecular clock of mitochondrial DNA. \emph{Journal of Molecular Evolution}
\bold{22}, 160--174.
}
\examples{
use_genetic_code("Vertebrate mitochondrial")
hky <- HasegawaKishinoYano(rep(0.25, 4), transition_rate = 2)
psi <- equilibrium_to_fitness(rep(0.05, 20), 1000)
m <- MutationSelection(1000, 1, hky, psi, scaling_type = "synonymous")
m
}
\seealso{\code{\link{equilibrium_to_fitness}}, \code{\link{simulate_over_phylogeny}}}
