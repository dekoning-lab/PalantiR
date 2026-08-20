# PalantiR 1.2.1

## Bug fixes

- `rescale_method = "exact"`: event positions through a stiff transient are
  now reported correctly. The exact rescaler's knot table marched at a fixed
  intrinsic-time step; when the entering forecast was far from the new mode's
  equilibrium the scaled class flux could start orders of magnitude above its
  equilibrium value (634x measured on a disjoint-profile switch), the whole
  boundary layer fell inside one knot interval, and the piecewise-linear
  inverse time change misplaced early events (4.6x too many class events
  reported by branch position 0.1). Event identities, total counts, and the
  exit forecast were unaffected. The knot grid now refines adaptively so that
  no knot interval delivers more than 1% of the branch's event budget, and
  the inverse time change locates intervals by binary search. Found by the
  engine validation notebook's transient test (`validation/`); guarded by
  regression test E1.

# PalantiR 1.2.0

Comprehensive fixes following a full source review
(31 verified findings; review notes ship with the analysis repository).

## Breaking changes

- `CoEvolution` models: the scaling constant is now computed over the actual
  codon-pair state space (previously it was meaningless, indexing only a 61x61
  corner of the pair matrix), so one unit of branch length is one expected
  substitution per codon; and instantaneous double substitutions have rate
  zero (previously they carried a dimensionally inconsistent product of two
  rates, making the normalised model depend on the absolute mutation rate).
  Time scales of all previous CoEvolution simulations were invalid.
- `MarkovModulatedMutationSelection`: the switching process is now a proper
  GTR, so the requested mode occupancy is realised (previously occupancy was
  always uniform and the requested weights instead inflated within-mode
  substitution rates by 1/w).
- `simulate_with_shared_time_heterogeneity`: two mode-assignment defects are
  fixed (an off-by-one that discarded the first sampled mode, and mode chains
  restarting from the parent's starting state); simulated mode paths differ
  from previous versions.
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

- Segment rescaler: multi-interval branches no longer over-run (previously up
  to ~3x re-simulation of a branch); per-segment memory reduced by ~50x;
  zero-length branches no longer crash (segments) or hang (exact).
- IntervalHistory segmentation: exact-multiple branch/segment pairs no longer
  drop the final segment or write out of bounds.
- Numerical: log-sum-exp offsets, fixation-probability guards for extreme
  selection, full-precision tree geometry in widgets.
