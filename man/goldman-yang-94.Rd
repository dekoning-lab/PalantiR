\name{goldman-yang-94}
\alias{GoldmanYang94}
\alias{GY94}
\alias{FEqual}
\alias{F1x4}
\alias{F3x4}
\alias{F61}
\alias{codon_frequencies}
\alias{GY94BranchModel}
\alias{GY94SiteModel}
\alias{DiscreteGY94}
\alias{simulate_gy94_site_model}
\alias{ZNYBranchSiteModel}
\alias{ZNYBranchSite}
\title{Goldman--Yang codon models and discrete site classes}
\usage{
FEqual()

F1x4(frequencies)

F3x4(frequencies)

F61(frequencies)

codon_frequencies(frequencies = NULL, frequency_model = NULL)

GoldmanYang94(omega = 1, kappa = 1, frequencies = NULL,
              frequency_model = NULL, scaling_type = "standard")

GY94(omega = 1, kappa = 1, frequencies = NULL,
     frequency_model = NULL, scaling_type = "standard")

GY94BranchModel(models, mode_phylogeny, start_mode = 0L)

GY94SiteModel(models, n_sites = NULL, site_classes = NULL,
              labels = names(models))

DiscreteGY94(models, n_sites = NULL, site_classes = NULL,
             labels = names(models))

simulate_gy94_site_model(phylogeny, site_model, sequence = NULL, rate = 1,
                         segment_length = 0.001, tolerance = 0.001,
                         rescale_method = "exact")

ZNYBranchSiteModel(mode_phylogeny, site_counts, omega0 = 0.2, omega2 = 2,
                   kappa = 2, frequencies = NULL, frequency_model = NULL,
                   foreground_frequencies = frequencies,
                   foreground_frequency_model = frequency_model,
                   scaling_type = "standard", start_mode = 0L)

ZNYBranchSite(mode_phylogeny, site_counts, omega0 = 0.2, omega2 = 2,
              kappa = 2, frequencies = NULL, frequency_model = NULL,
              foreground_frequencies = frequencies,
              foreground_frequency_model = frequency_model,
              scaling_type = "standard", start_mode = 0L)
}
\arguments{
\item{omega}{Positive non-synonymous/synonymous rate ratio.}
\item{kappa}{Positive transition/transversion rate ratio.}
\item{frequencies}{A frequency specification. For \code{GoldmanYang94} and
\code{codon_frequencies}, this may be a \code{CodonFrequencies} object, four
nucleotide frequencies, a position-specific 3 by 4 or 4 by 3 matrix, a
length-12 vector, or one value per sense codon. For \code{F1x4}, \code{F3x4}
and \code{F61}, it is the corresponding nucleotide or codon representation.}
\item{frequency_model}{One of \code{"FEqual"}, \code{"F1x4"},
\code{"F3x4"}, or \code{"F61"}. Case and punctuation are ignored. If omitted,
the representation is inferred from the dimensions of \code{frequencies};
\code{NULL} selects FEqual.}
\item{scaling_type}{Branch-length normalization. \code{"standard"} is an alias
of \code{"substitution"}. \code{"synonymous"} is the dS gauge;
\code{"dS"} and \code{"ds"} are explicit aliases. The other scaled choices
are \code{"synonymous-per-codon"} and \code{"non-synonymous"}.
\code{"none"} leaves the generator unscaled.}
\item{models}{For \code{GY94BranchModel}, a list of GY94 models indexed by the
mode labels. For \code{GY94SiteModel}, a list whose elements are either GY94
models or \code{GY94BranchModel} specifications.}
\item{mode_phylogeny}{A \code{Phylogeny} of type \code{"mode"}, with the same
topology as the branch-length phylogeny. Modes are zero-based indices into a
branch model's \code{models}.}
\item{start_mode}{Zero-based index of the model in force at the root.}
\item{n_sites}{One positive whole-number site count per class. Supply exactly
one of \code{n_sites} and \code{site_classes}.}
\item{site_classes}{One class assignment per output site. Character or factor
values refer to \code{labels}; numeric values are one-based indices into
\code{models}. Every model must be assigned at least one site.}
\item{labels}{One unique, non-empty character label per site class.}
\item{phylogeny}{A branch-length \code{Phylogeny} over which to simulate.}
\item{site_model}{A \code{GY94SiteModel} from \code{GY94SiteModel} or
\code{ZNYBranchSiteModel}.}
\item{sequence}{Optional root codon \code{Sequence} with one state per output
site. If omitted, each site's root is sampled from the equilibrium of its
class's \code{start_mode} model.}
\item{rate}{Scalar multiplier applied to every branch length.}
\item{segment_length, tolerance, rescale_method}{Transient-rescaling controls
passed to \code{simulate_over_interval_phylogeny} for branch-heterogeneous
classes. See \code{\link{simulate_over_interval_phylogeny}}.}
\item{site_counts}{Four positive whole-number counts, in class order 0, 1, 2a,
and 2b.}
\item{omega0}{Background purifying-selection omega used by classes 0 and 2a.}
\item{omega2}{Foreground omega used by classes 2a and 2b.}
\item{foreground_frequencies, foreground_frequency_model}{Optional frequency
specification for foreground branches. By default the background specification
is reused.}
}
\description{
Construct and simulate Goldman--Yang 1994 codon models. The interface supports
standard codon-frequency parameterizations, branch-specific model changes,
arbitrary discrete site classes, and the four-class Zhang--Nielsen--Yang
branch-site model.
}
\details{
\strong{Generator.} For distinct sense codons \eqn{i} and \eqn{j}, the
unscaled rate is zero unless they differ at exactly one nucleotide. Otherwise
it is \eqn{\pi_j} multiplied by \eqn{\kappa} for a transition, by \eqn{\omega}
for a non-synonymous change, or by both factors when both conditions apply.
Diagonal entries are negative row sums. Both \code{omega} and \code{kappa} must
be finite and strictly positive.

\strong{Frequency adapters.} \code{FEqual} assigns equal mass to the sense
codons of the active genetic code. \code{F1x4} multiplies one set of nucleotide
frequencies across the three positions. \code{F3x4} uses a different set at
each codon position. Unnamed nucleotide input uses PalantiR's order T, C, A, G;
named input is reordered. \code{F61} takes one positive value per sense codon
(the actual count follows the active genetic code) and should preferably be
named by codon. Stop codons are excluded and all four adapters normalize over
the remaining sense codons. A \code{CodonFrequencies} object is tied to the
genetic code under which it was constructed and cannot be reused after the
active code changes.

\strong{Scaling.} For \code{"synonymous"}, PalantiR constructs the neutral
\eqn{\omega=1} generator \eqn{\widetilde Q^{(1)}} and divides the selected
generator by
\eqn{S_{dS}=-(1/3)\sum_i\pi_i\widetilde q^{(1)}_{ii}}. The neutral total rate
is then three substitutions per codon, and branch length is dS per synonymous
opportunity. This denominator is independent of the simulated \eqn{\omega}.
The aliases \code{"dS"} and \code{"ds"} select the same gauge, stored as
\code{"dS"}.

For event-count scaling, let \eqn{\mathcal A} be the set of ordered codon transitions
that determines scaling and let \eqn{\widetilde Q} denote the unscaled
generator. Its stationary rate is
\eqn{S_{\mathcal A}=\sum_i\pi_i\sum_{j:(i,j)\in\mathcal A}
\widetilde q_{ij}}, where \eqn{\mathcal A} contains all,
synonymous-per-codon, or non-synonymous substitutions. A standalone model is divided by
\eqn{S_{\mathcal A}}. \code{"standard"} is accepted as an alias and stored as
the canonical value \code{"substitution"}. The explicit
\code{"synonymous-per-codon"} name preserves PalantiR's earlier
\code{"synonymous"} behavior.

For a discrete mixture with class proportions \eqn{p_k}, every class shares
the mixture-wide denominator
\eqn{\overline S_{\mathcal A}=\sum_k p_k S_{\mathcal A,k}}. Class \eqn{k}
therefore has stationary scaled rate
\eqn{S_{\mathcal A,k}/\overline S_{\mathcal A}}, and their weighted mean is
one. For a branch-site model, a separate mixture-wide denominator is computed
for each branch type, and all classes active on that branch type share it. This
retains relative rates among site classes and follows the PAML convention. All
classes and modes in one mixture must use the same scaling type.

\strong{Branch heterogeneity.} A \code{GY94BranchModel} may change omega, kappa,
or stationary codon frequencies among branches. Following a frequency change,
the child process starts from the distribution inherited from its parent; it is
not redrawn from the new equilibrium. Event-count gauges use the normal
transient-rescaling machinery to preserve the selected event budget during
this approach. The dS gauge uses the same machinery with the neutral-reference
total outflux per codon divided by three as its clock functional. This
preserves neutral nucleotide time without imposing a realized synonymous-event
count.

\strong{Discrete classes.} \code{GY94SiteModel} is the general representation.
Each class may contain one homogeneous GY94 model or one branch-heterogeneous
\code{GY94BranchModel}. Count-based assignments place classes in list order;
explicit assignments are returned in the requested site order. All
branch-heterogeneous classes in one mixture must use the same mode phylogeny and
number of branch models so that one common denominator can be applied on each
branch type.

\strong{Zhang--Nielsen--Yang classes.} \code{ZNYBranchSiteModel} constructs
class 0 as \eqn{\omega_0 \rightarrow \omega_0}, class 1 as
\eqn{1 \rightarrow 1}, class 2a as \eqn{\omega_0 \rightarrow \omega_2}, and
class 2b as \eqn{1 \rightarrow \omega_2}. The mode tree should label background
branches 0 and foreground branches 1.
}
\value{
\code{FEqual}, \code{F1x4}, \code{F3x4}, \code{F61}, and
\code{codon_frequencies} return a named numeric vector of class
\code{CodonFrequencies}, normalized over the active code's sense codons.

\code{GoldmanYang94} and \code{GY94} return a
\code{c("GoldmanYang94", "SubstitutionModel")} object containing the
stationary distribution, normalized transition matrix, simulation sampling
matrix, omega, kappa, frequency representation, scaling constant, scaling type,
and genetic code. dS-scaled models also report the raw neutral total and
synonymous rates and the Q-weighted synonymous opportunities used to define
the gauge.

\code{GY94BranchModel} returns a branch-process specification and
\code{GY94SiteModel}, \code{DiscreteGY94}, \code{ZNYBranchSiteModel}, and
\code{ZNYBranchSite} return a \code{GY94SiteModel}. The latter records class
weights, the common scaling denominator for every branch mode, and each
class's stationary rate relative to that denominator.

\code{simulate_gy94_site_model} returns a \code{Simulation}. Its
\code{site_classes} data frame contains every zero-based site, its character
\code{site_class}, and its zero-based \code{class_index}, including sites with
no events. The same table is attached to the alignment. Each row of
\code{substitutions} contains \code{site_class} and \code{class_index};
\code{site_model} and \code{site_class_models} retain the generating
specification.
}
\references{
Goldman, N. and Yang, Z. (1994). A codon-based model of nucleotide substitution
for protein-coding DNA sequences. \emph{Molecular Biology and Evolution}
\bold{11}, 725--736.

Yang, Z., Nielsen, R., Goldman, N. and Pedersen, A.-M. K. (2000).
Codon-substitution models for heterogeneous selection pressure at amino acid
sites. \emph{Genetics} \bold{155}, 431--449.

Yang, Z. and Nielsen, R. (2002). Codon-substitution models for detecting
molecular adaptation at individual sites along specific lineages.
\emph{Molecular Biology and Evolution} \bold{19}, 908--917.

Zhang, J., Nielsen, R. and Yang, Z. (2005). Evaluation of an improved branch-
site likelihood method for detecting positive selection at the molecular
level. \emph{Molecular Biology and Evolution} \bold{22}, 2472--2479.
}
\examples{
use_genetic_code("Standard nuclear")

pi <- F1x4(c(T = 0.3, C = 0.2, A = 0.25, G = 0.25))
purifying <- GY94(omega = 0.2, kappa = 2, frequencies = pi)
neutral <- GY94(omega = 1, kappa = 2, frequencies = pi)

classes <- GY94SiteModel(
    models = list(purifying = purifying, neutral = neutral),
    n_sites = c(4, 2))

tree <- Phylogeny(system.file("extdata", "one_branch.newick",
                              package = "PalantiR"))
set_palantir_seed(1)
sim <- simulate_gy94_site_model(tree, classes)
sim$site_classes
}
\seealso{
\code{\link{simulate_over_phylogeny}},
\code{\link{simulate_over_interval_phylogeny}},
\code{\link{Phylogeny}},
\code{\link{GeneticCode}}
}
