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
