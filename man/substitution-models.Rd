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
                  scaling_type = "synonymous", fixation = "diffusion")

CoEvolution(population_size, mutation_rate, nucleotide_model, fitness_1,
            fitness_2, delta, scaling_type = "synonymous",
            fixation = "diffusion")

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
\item{fixation}{The fixation probability of a new mutant: \code{"diffusion"}
(Kimura's eq. 10, the default) or \code{"weak-selection"} (his eq. 11, as in
Halpern and Bruno, PhyloBayes MutSel and swMutSel). See Details.}
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
Mutation-selection and co-evolution models also store \code{population_size}
and \code{fixation}.
}
\description{
Construct substitution models. \code{MutationSelection} implements the
mutation-selection formulation of Halpern and Bruno (1998) over codons, in which
amino-acid fitnesses and an effective population size determine fixation
probabilities and hence substitution rates.
}
\details{
\strong{Fixation probability.} A single-nucleotide move from codon \eqn{i} to
codon \eqn{j} has rate \eqn{2N\mu_{ij}u(s)}, where \eqn{N} is the number of
diploid individuals, \eqn{\mu_{ij}} the mutation rate and
\eqn{s = f_j - f_i} the difference between the fitnesses of the two amino
acids, which is the homozygote selection coefficient (Kimura's \eqn{2s}).
\code{"diffusion"} uses Kimura's (1962) diffusion result,
\deqn{u(s) = \frac{1 - e^{-s}}{1 - e^{-2Ns}},}{u(s) = (1 - exp(-s)) / (1 - exp(-2Ns)),}
and \code{"weak-selection"} its weak-selection limit,
\deqn{u(s) = \frac{s}{1 - e^{-2Ns}}.}{u(s) = s / (1 - exp(-2Ns)).}
Under weak selection the rate is \eqn{\mu_{ij} S / (1 - e^{-S})} with
\eqn{S = 2Ns}, a function of the population-scaled coefficient alone. This is
the substitution model of Halpern and Bruno (1998), PhyloBayes MutSel
(Rodrigue et al. 2010) and swMutSel (Tamuri et al. 2012); use it to simulate
data for those programs, whose fitted profiles are scaled fitnesses
\eqn{2Nf}. Under the diffusion form the rate also depends on \eqn{N}: the two
rates differ by a factor of about \eqn{1 - S/4N}, and the diffusion form is
the closer approximation to Wright-Fisher fixation at small \eqn{N}.

The stationary distribution of each codon is proportional to the mutational
weight of its nucleotides times \eqn{e^{Ef}}, where detailed balance gives
\eqn{E = 2N - 1} under the diffusion form and \eqn{E = 2N} under weak
selection. \code{equilibrium} is exactly this distribution. Before version
1.4.0 it used \eqn{E = 2N} under the diffusion form as well, and so was not
the stationary distribution of \code{transition}; each frequency ratio was off
by a factor \eqn{e^{f_j - f_i}}.

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
Kimura, M. (1962) On the probability of fixation of mutant genes in a
population. \emph{Genetics} \bold{47}, 713--719.

Halpern, A. L. and Bruno, W. J. (1998) Evolutionary distances for protein-coding
sequences: modeling site-specific residue frequencies.
\emph{Molecular Biology and Evolution} \bold{15}, 910--917.

Hasegawa, M., Kishino, H. and Yano, T. (1985) Dating of the human-ape splitting
by a molecular clock of mitochondrial DNA. \emph{Journal of Molecular Evolution}
\bold{22}, 160--174.

Rodrigue, N., Philippe, H. and Lartillot, N. (2010) Mutation-selection models
of coding sequence evolution with site-heterogeneous amino acid fitness
profiles. \emph{Proceedings of the National Academy of Sciences USA}
\bold{107}, 4629--4634.

Tamuri, A. U., dos Reis, M. and Goldstein, R. A. (2012) Estimating the
distribution of selection coefficients from phylogenetic data using sitewise
mutation-selection models. \emph{Genetics} \bold{190}, 1101--1115.
}
\examples{
use_genetic_code("Vertebrate mitochondrial")
hky <- HasegawaKishinoYano(rep(0.25, 4), transition_rate = 2)
psi <- equilibrium_to_fitness(rep(0.05, 20), 1000)
m <- MutationSelection(1000, 1, hky, psi, scaling_type = "synonymous")
m
}
\seealso{\code{\link{equilibrium_to_fitness}}, \code{\link{simulate_over_phylogeny}}}
