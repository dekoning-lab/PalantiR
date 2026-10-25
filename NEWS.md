# PalantiR 1.4.0

## Stationary distribution of the mutation-selection generator

- `MutationSelection()` and `CoEvolution()` simulate with Kimura's (1962)
  diffusion fixation probability, u(s) = (1 - e^-s) / (1 - e^-2Ns) (his
  eq. 10, with s the homozygote selection coefficient). Detailed balance for
  this generator gives a stationary distribution proportional to the
  mutational weight times exp((2N - 1) f). The reported `equilibrium`, and
  hence the root sequences sampled from it, used exp(2N f), the stationary
  distribution of the weak-selection form (eq. 11) instead. Each frequency
  ratio was off by a factor exp(f_j - f_i), which is exp(S/2N) for a scaled
  difference S: with scaled fitnesses spanning 20 units, amino-acid
  frequencies were off by up to 9% at N = 100, 1% at N = 1000 and 1e-5 at
  N = 10^6. Descendant sequences relaxed from the reported distribution
  towards the true one along the tree. `equilibrium` now uses 2N - 1 and is
  the stationary distribution of `transition` to rounding error.
- `equilibrium_to_fitness()` now divides by the same exponent, 2N - 1, so
  that a model built on its output has the requested stationary distribution.
  Its default output changes: every fitness moves away from the fittest
  amino acid's by the factor 2N / (2N - 1). Simulations at large N are
  essentially unaffected (a 5e-7 relative change in every selection
  coefficient at N = 10^6).

## New features

- `MutationSelection()`, `CoEvolution()` and `equilibrium_to_fitness()` take
  `fixation = "diffusion"` (the default, as above) or
  `fixation = "weak-selection"`, Kimura's eq. (11),
  u(s) = s / (1 - e^-2Ns). Its substitution rate, mu S / (1 - e^-S) with
  S = 2Ns, depends on the population-scaled coefficient alone and is the
  model of Halpern and Bruno (1998), PhyloBayes MutSel and swMutSel; its
  stationary distribution carries exp(2N f). Use it to simulate data for
  those programs exactly. `equilibrium_to_fitness(..., fixation =
  "weak-selection")` reproduces the pre-1.4.0 default output bit for bit.
  Models store the choice in `fixation`.

# PalantiR 1.3.0

## Corrected C50 profiles

- `psi_c50_1` had its V and W columns exchanged: the column documented as V
  held the W frequencies and vice versa, so profile 2 (41% I, 43% V in the
  source) read as 41% I, 43% W. The swap was present when the file was added
  to Palantir-Perl in 2015 and was copied into PalantiR unchanged. The data are
  now relabelled, every row reproduces Table 4 of the supplementary material
  of Le, Gascuel & Lartillot (2008, *Bioinformatics* 24:2317-2323), and the
  matrix has amino-acid column names. The raw `inst/extdata/psi-c50-1.txt` is
  corrected the same way. Simulations that used `psi_c50_1` in earlier
  versions used profiles with V and W exchanged; rerun any whose conclusions
  depend on the profiles' composition.
- New data set `psi_c50_2`: the second independent C50 estimate (Table 5 of the
  same supplement), with raw file `inst/extdata/psi-c50-2.txt`.
- New data set `psi_c50_iqtree`, with weights `psi_c50_iqtree_weights`: the C50
  profiles and mixture weights built into IQ-TREE 2, reordered to PalantiR's
  amino-acid order. These are a different estimate of C50 from Tables 4 and 5
  (no profile matches a table row after rounding) and are not interchangeable
  with `psi_c50_1` or `psi_c50_2`.
- `?datasets` now cites the source tables and describes the raw files'
  trailing 21st field correctly (a zero for Palantir-Perl's stop symbol, not a
  trailing tab).

## Corrected branch-length semantics

- `scaling_type = "synonymous"` now implements dS per synonymous opportunity,
  using a neutral total rate of three substitutions per codon. `"dS"` and
  `"ds"` are explicit aliases; models store the canonical value `"dS"`.
  The previous one-synonymous-event-per-codon convention remains available as
  `"synonymous-per-codon"`.
- GY94 dS scaling is computed from the `omega = 1` reference and is independent
  of the simulated omega. Mutation--selection dS scaling is computed from the
  mutation-only neutral reference, so site-specific fitness profiles sharing a
  mutation model retain the same clock.
- Time-heterogeneous dS simulations use the same exact or segmented transient
  rescaler as other gauges, with a dS-specific clock functional: the
  neutral-reference total outflux per codon divided by three. Fitness-profile
  and population-size shifts therefore preserve the intended dS budget without
  imposing one realized synonymous event per codon.

## New features

- `GoldmanYang94()` (also `GY94()`) implements the Goldman--Yang 1994 codon
  model with `omega`, `kappa`, and Fequal, F1x4, F3x4, or F61 stationary-codon
  frequencies. It supports standard/all-substitution, synonymous, and
  non-synonymous branch-length scaling; `"standard"` is accepted as an alias
  and stored canonically as `"substitution"`.
- `GY94BranchModel()` assigns different GY94 processes to branches using a mode
  tree, including changes in omega, kappa, or stationary codon frequencies.
  The existing transient rescalers preserve the requested branch-length
  currency when the process changes.
- `GY94SiteModel()` and `simulate_gy94_site_model()` provide arbitrary discrete
  site classes, each of which may be time homogeneous or branch heterogeneous.
  Classes share one mixture-weighted scaling denominator, with a separate
  denominator for each branch type in branch-site models, following the PAML
  convention. Relative rates among classes are retained; individual component
  rates need not equal one.
  Site-class assignments are retained in the simulation, alignment metadata,
  substitution table, and interactive alignment viewer.
- `ZNYBranchSiteModel()` constructs the four standard Zhang--Nielsen--Yang
  branch-site classes (0, 1, 2a, and 2b) as a special case of the general
  discrete-class simulator.
- A rerunnable R-kernel notebook now performs the complete 50-taxon,
  5,000-codon branch-heterogeneous GY94 simulation directly with PalantiR under
  synonymous scaling. It explores the live alignment and history widgets,
  reconstructs Goldman--Yang dN/dS from all 762,540 recorded events, visualizes
  a complete site history, and compares realized with exact rates. The original
  Python notebook remains as an independent audit of the saved production
  artifact.

# PalantiR 1.2.2

- `rescale_method = "exact"` is now the default for
  `simulate_over_interval_phylogeny`. The exact time change replaces the
  branch-segmentation scheme it approximated: it delivers the branch's event
  budget identically rather than to within the segment tolerance, runs
  several-fold faster, and (since 1.2.1) resolves stiff transients correctly.
  The segment scheme remains available with `rescale_method = "segments"`.

# PalantiR 1.2.1

## Bug fixes

- `rescale_method = "exact"`: event positions through a stiff transient are
  now reported correctly. Event identities, total counts, and the
  exit forecast of the transient probability distribution were unaffected.
  The knot grid now refines adaptively so that no knot interval delivers more
  than 1% of the branch's event budget, and the inverse time change locates
  intervals by binary search.

# PalantiR 1.2.0

Comprehensive fixes (31 verified findings; see review notes).

## Breaking changes

- `CoEvolution` models: the scaling constant is now computed over the actual
  codon-pair state space (previously it incorrectly indexed only a corner of
  the pair matrix), so one unit of branch length is one expected
  substitution per codon; and instantaneous double substitutions have rate
  zero. Time scales of all previous CoEvolution simulations are suspect and
  should be rerun.
- `MarkovModulatedMutationSelection`: the switching process now follows a
  GTR model.
- `simulate_with_shared_time_heterogeneity`: two mode-assignment defects are
  fixed (an off-by-one error previoulsy discarded the first sampled mode);
  simulated mode paths will differ from previous versions.
- Guide/mode trees must now be constructed with `Phylogeny(..., type = "mode")`
  and are validated; passing a branch-length tree where a mode tree is
  expected is an error instead of a silent single-mode simulation.
- Input validation throughout: fitness and delta dimensions, positive
  population sizes and mutation rates, Newick sanity, genetic code names
  (`use_genetic_code` now validates against `genetic_code_names()`).
- Datasets are now proper `.rda` objects (`data(amino_acid_fitness_N_1000)`
  et al. previously loaded as unusable single-column tables); `psi-c50-1` is
  renamed `psi_c50_1`.
- `sim$intervals` is `NULL` (not `0`) for time-homogeneous simulations.
- `join()` offsets each simulation's 0-based `site` column so joined
  substitution tables no longer collide.

## New features

- `set_palantir_seed(seed)`: the simulation engine is now reproducible. Note
  that base R `set.seed()` does not affect the engine.
- `genetic_code_names()`: the valid names for `use_genetic_code()`.
- `equilibrium_to_fitness(pi, N, nucleotide_equilibrium = ...)`: optional
  mutational-weight correction so the realised amino-acid marginal equals the
  target; without the argument, behaviour is unchanged. `N` may now be any
  positive real number (it was previously truncated to an integer silently).
- `plot(phylogeny, mode_phylogeny)` works positionally, and widget options
  passed through `plot()` now reach the widget.

## Fixes

- Segment rescaler: multi-interval branches no longer over-run or require
  repeat re-simulation; per-segment memory reduced by ~50x; zero-length
  branches no longer crash or hang.
- IntervalHistory segmentation: exact-multiple branch/segment pairs no longer
  drop the final segment or write out of bounds.
- Numerical: log-sum-exp offsets, fixation-probability guards for extreme
  selection, full-precision tree geometry in widgets.
