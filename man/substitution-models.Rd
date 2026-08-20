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
\code{"synonymous"}, \code{"non-synonymous"}, \code{"substitution"} or
\code{"none"}. See Details.}
\item{switching_model}{An exchangeable model governing switching between modes.}
\item{mutation_selection_models}{List of mutation-selection models, one per mode.}
}
\value{
A list of class \code{SubstitutionModel} with elements \code{equilibrium},
\code{transition} (the normalised rate matrix), \code{sampling} (its cumulative
form, used by the simulators), \code{scaling} (the normalising constant applied),
\code{scaling_type}, \code{n_states} and \code{type}.
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

\strong{Scaling.} The rate matrix is divided by a normalising constant so that
one unit of branch length corresponds to one expected substitution of the chosen
class per codon. \code{"synonymous"} makes branch lengths count synonymous
substitutions, \code{"non-synonymous"} amino-acid-changing substitutions,
\code{"substitution"} all substitutions, and \code{"none"} applies no scaling.
Choose the class matching the units in which the guide tree's branch lengths were
estimated. An unrecognised name raises an error; note that the option called
\code{"standard"} in the 2016 version of this package is \code{"substitution"}
here.
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
