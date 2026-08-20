# HANDOFF

Notes from one agent to another about changes that need work in files the
writing agent does not own.

**Outstanding edit needed in `NAMESPACE`:** add `S3method("[", Alignment)`.
See "From the R/widget agent" section 1 below — `tests_review/MIN9.R` fails
without it.

## From the C++/API agent (C7, C8, M2, M3, M4, M5, M7 C++ half, MIN6, MIN7)

Files changed by this agent: `src/Models.cpp`, `src/Palantir.cpp`,
`src/Sequence.cpp`, `src/RcppPalantir.cpp`, `src/RcppPalantir.hpp`,
`src/RcppExports.cpp`, `src/Palantir_Core/Palantir.{cpp,hpp}`,
`src/Palantir_Core/Util.hpp`, `R/RcppExports.R`.

### 1. BREAKING: guide trees must now be built with `type = "mode"` (M5)

`simulate_over_interval_phylogeny()` and `phylogeny_to_intervals()` now reject a
guide tree whose `type` is not `"mode"`, and all four simulators reject a
primary tree whose `type` is not `"phylogeny"`. `Phylogeny()` also rejects a
`type` other than those two.

Known call sites outside the package that this breaks (not mine to edit):

- `repro/simulate/genSimNonSyn.R:86` — `md <- Phylogeny(guideTree)` must become
  `md <- Phylogeny(guideTree, type = "mode")`.
- `analysis/genSimHold.R`, `analysis/engine_tests/shim.R`,
  `analysis/engine_tests/t6_rescale2016.R` call a different, older
  `Phylogeny(path, attribute_file =, attr =)` signature and do not run against
  this package version at all; no action implied.

`R/Phylogeny.R`'s `plot.Phylogeny` already requires `mode_phylogeny$type ==
"mode"`, so it is unaffected. Any README/vignette/man example that builds a
guide tree needs the argument added.

### 2. `genetic_code_names()` is available for R-side validation (M7)

New exported C++ function, no arguments, returns a character vector of the 18
valid genetic code names (read straight from the C++ `GeneticCode_table`, so it
cannot drift). `use_genetic_code()` should validate against it, e.g.

```r
use_genetic_code <- function(name) {
    valid <- genetic_code_names()
    if (!is.character(name) || length(name) != 1L || is.na(name) || !(name %in% valid)) {
        stop("Unknown genetic code ...; valid names are: ", paste(sort(valid), collapse = ", "))
    }
    .globals$genetic_code_name <- name
}
```

`man/genetic-code.Rd` also has an invalid example name (`"Standard"`; the real
name is `"Standard nuclear"`) and should document `genetic_code_names()`.

### 3. `sim$intervals` is now really `NULL` for homogeneous simulations (MIN6)

`simulate_over_phylogeny()` used to return the integer `0` there (Rcpp wrapping
the C++ `NULL` pointer literal). It is now `R_NilValue`. Anything that branched
on `is.null(sim$intervals)` — `plot.Simulation` passes it straight to
`PhyloPlot()`, and `join()` collects it — will now take the other branch. Worth
a look from whoever owns `R/Palantir.R`, `R/Simulation.R` and the widget JS
(this interacts with M9).

### 4. New/changed exported signatures needing man pages (`man/` is not mine)

- `set_palantir_seed(seed)` — NEW. Seeds the single engine `mt19937`. Document
  prominently that **R's `set.seed()` does not affect PalantiR's engine**; this
  is the only way to make a simulation reproducible. `seed` must be a whole
  number in `[0, 4294967295]`.
- `genetic_code_names()` — NEW, see above.
- `equilibrium_to_fitness(equilibrium, population_size, nucleotide_equilibrium = NULL)`
  — third argument is NEW and optional. Omitted, the function is bit-identical
  to before (the degeneracy-tilted inversion; `man/equilibrium_to_fitness.Rd`
  should keep its current description for that case and gain the caveat).
  Supplied (length 4, the nucleotide equilibrium of the mutation model), the
  mutational weight of each amino acid's codons is divided out so that the
  amino-acid marginal of the resulting `MutationSelection` equilibrium equals
  `equilibrium` exactly. Verified round trip: TV distance ~2e-14.
  `\value` should also say it returns a 20x1 matrix, not a vector.
- `MutationSelection`, `CoEvolution`, `equilibrium_to_fitness`:
  `population_size` must now be a whole number >= 1; `mutation_rate` must be
  > 0; `sample_sequence`'s `length` must be a whole number >= 1;
  `fitness`/`fitness_1`/`fitness_2` must be length 20 and `delta` must be 20x20.
  `man/substitution-models.Rd` and `man/sequence.Rd` should state these.

### 5. Not done here, still open

- The 2-argument `equilibrium_to_fitness` tilt is a *manuscript caveat*, not a
  bug that was silently repaired: the published pipeline's numbers are
  unchanged by these edits.
- `man/simulate.Rd` still omits `rescale_method`, and the two shared-*
  simulators are undocumented (from the minor list; not touched).

## From the simulation-engine agent (C1, C2, C3, C4a, C4b, M11, M12, MIN1-MIN4)

Files changed by this agent: `src/Palantir_Core/Simulate.cpp`,
`src/Palantir_Core/IntervalHistory.cpp`, `src/Palantir_Core/IntervalHistory.hpp`.
No exported signature changed.

### 1. NEW ERROR: `rescale_method = "exact"` requires a scaled model class (MIN1)

`sequence_over_intervals` now throws when `rescale_method` is `"exact"` and the
models were built with `scaling_type = "none"`, instead of silently running the
segment scheme instead. It throws up front, whether or not the rescaler would
have engaged. This reaches R through `simulate_over_interval_phylogeny()`.

`man/simulate.Rd` still does not document `rescale_method` at all (already on
the minor list); when it is written up it should say that `"exact"` needs
models built with a `scaling_type` other than `"none"`.

### 2. The mode chain from `simulate_with_shared_time_heterogeneity` has changed (C2, C3)

`switching_poisson` was giving each interval the mode entered at the *next*
switch, and never wrote a branch's end-of-branch mode back, so every child
restarted from the root's `start_mode`. Both are fixed. `sim$intervals` from
that simulator therefore reports a different — and now Markov — mode chain.
Nothing about the shape of the output changed, but anything that consumed those
per-site interval tables (the widget path in M9, `join()`) is now looking at
different numbers. Consider it a behaviour change worth a NEWS entry.

### 3. The exact rescaler is unreachable from the shared simulators

`simulate_with_shared_substitution_heterogeneity()` and
`simulate_with_shared_time_heterogeneity()` in `src/Palantir.cpp` (not mine)
do not expose `rescale_method` and so always use the segment scheme, at their
default `segment_length = 1e-4`. That is the configuration M12 was about. The
refactor cuts its memory cost by roughly two orders of magnitude, but the exact
method is still both faster and cheaper, and `sequence_over_intervals` already
takes the argument. Whoever owns `src/Palantir.cpp` may want to add a
`rescale_method` parameter to both, defaulting to `"segments"` so nothing
changes for existing callers.

### 4. Zero-length branches now simulate nothing rather than crashing (C4a, C4b)

An interval with `finish <= start` is skipped by the rescaler: no segment is
pushed, and the child inherits the forecast and generator that entered the
branch. Such a branch carries no substitutions (it never could have; it used to
segfault in segment mode and run forever in exact mode). `sim$intervals` still
lists the zero-length interval, because that table comes from `to_intervals` /
`switching_*`, not from the rescaler.

### 5. `tests_review/MIN2.R`, `MIN3.R` and `MIN4.R` assert on source text

Those three fixes are, respectively, below Monte-Carlo resolution, a removal of
unreachable code, and a change to C++ argument *names* — none of them can be
distinguished from outside the compiled library, so the tests assert on the
contents of `src/Palantir_Core/Simulate.cpp` and
`src/Palantir_Core/IntervalHistory.{cpp,hpp}` (plus a behavioural sanity check
each). If anyone reformats or moves those files, update the markers those tests
grep for: `// done rescaling - rest of branch`, `else { // entire interval`,
`j + 1 == K`, and the two-`const double` constructor signature.

## From the docs/data agent (D1-D4)

Files changed by this agent: `man/*.Rd`, `README.md`, `DESCRIPTION`,
`data/`, `inst/extdata/`, `docs/*.Rmd`.

### 1. The datasets are now .rda, and one was renamed (D1)

`data/` held `amino_acid_equilibrium.csv`, `amino_acid_fitness_N_1000.csv` and
`psi-c50-1.txt`; `data()` read them with `read.table`'s `sep = ";"` convention
and delivered unusable single-column objects. Those three raw files have moved
to `inst/extdata/` (unchanged, byte for byte) and `data/` now holds
`amino_acid_equilibrium.rda`, `amino_acid_fitness_N_1000.rda` and
`psi_c50_1.rda`. `LazyData: true` is unchanged.

- `amino_acid_equilibrium`, `amino_acid_fitness_N_1000`: 6 x 20 numeric
  matrices, columns named for the amino acids.
- `psi_c50_1`: 50 x 20 numeric matrix, no dimnames. **Renamed** from the
  non-syntactic `psi-c50-1`; the raw file keeps its old name.
  Every line of the raw file ends in a tab, which `read.delim` turns into a
  21st, all-zero column; the packaged matrix has it removed.

Consequences outside this agent's files: anything reading `data/*.csv` by path
must now read `inst/extdata/` (or use `system.file`/`data()`). The `docs/*.Rmd`
vignettes were repointed at `../inst/extdata/`; their committed `.html`
renderings are stale on this point and were not regenerated.

### 2. Man pages now promise behaviour that lives in R/ (D2)

`man/genetic-code.Rd` documents `genetic_code_names()` and states that
`use_genetic_code()` **rejects a name that is not in that list**. The C++ half
(`genetic_code_names()`) exists; the R-side check in `R/Palantir.R` is not this
agent's to write (see section 2 of the C++/API agent's notes above). If that
validation does not land, the man page overstates what the function does.

`man/widgets.Rd` documents `sites` on the widgets, and `man/simulate.Rd` the
`site` column of the substitution table, as 0-based -- describing current
behaviour, not requesting a change.

## From the R/widget agent (M7 R half, M8, M9, M10, MIN8-MIN12)

Files changed by this agent: `R/Phylogeny.R`, `R/Palantir.R`, `R/Simulation.R`,
`R/Alignment.R`, `inst/htmlwidgets/AlignmentPlot.js`,
`inst/htmlwidgets/PhyloPlot.js`, `inst/htmlwidgets/PhyloPlot/Tooltip.js`.
`R/SubstitutionModel.R` needed nothing.

### 1. BLOCKING: `NAMESPACE` needs one more S3method line (MIN9)

`R/Alignment.R` now defines

```r
`[.Alignment` <- function(x, i, j, ..., drop = FALSE)
```

so that `alignment[, 1:3]` keeps its class, its `type` attribute and its taxon
rownames instead of decaying to a plain character matrix (which `as.fasta()` and
`plot()` then refused to dispatch on). The name starts with `[`, so
`exportPattern("^[[:alpha:]]+")` does NOT match it and internal `[` dispatch from
the user's workspace will not find it. **Add to `NAMESPACE`, next to the other
S3method lines:**

```r
S3method("[", Alignment)
```

Without that line `tests_review/MIN9.R` fails with exactly the message the
pre-fix library gave ("no applicable method"), and the fix is inert for users.
Everything else in this agent's batch is independent of it.

Two internal call sites were adjusted for the new method, because a row slice of
an Alignment is now a one-row Alignment rather than a bare character vector:
`.normalize()` returns `unclass(alignment)` for the codon/nucleotide case (the
other branches already returned plain matrices from `apply()`), and
`as.fasta.Alignment()` works on `unclass(x)`. Without those two, `AlignmentPlot()`
died in `toJSON()` with "No method asJSON S3 class: Alignment".

### 2. `plot.Phylogeny()` changed signature (M8)

`plot.Phylogeny(x, ..., mode_phylogeny = NULL)` is now
`plot.Phylogeny(x, mode_phylogeny = NULL, ...)`, and `...` is forwarded to
`PhyloPlot()` in both branches. The documented positional call
`plot(phylogeny, mode_phylogeny)` now works (it silently rendered a plain tree
before), and widget options passed to `plot()` now take effect. `man/` should
document the second argument and note that any `PhyloPlot()` option can be passed
through. Nothing that called it by name changes.

### 3. `join()` output changed shape (M10)

- `substitutions$site` of the *k*th simulation is now offset by the number of
  alignment columns contributed by simulations 1..*k*-1, so `site` indexes the
  joined alignment. Anything downstream that assumed the old (colliding) indices
  was wrong before.
- `intervals` is now one interval table **per site of the joined alignment**,
  named by the 0-based site index (`joined$intervals[["3"]]` belongs to alignment
  column 4), or `NULL` if no input simulation had intervals. It used to be one
  entry per *simulation*, which nothing could consume.
- `models` was `.get_member(sims, "model")`, which is `NULL` for every
  heterogeneous simulation (they store `models`, not `model`), so joining two
  interval simulations silently produced `models = list(NULL, NULL)`. It now
  takes `model` or `models`, one entry per simulation. This is not on the
  findings list; it is a one-line consequence of touching the same list.

`man/` for `join()` (if any) should state that the joined `site` column indexes
the joined alignment and describe the per-site `intervals` list.

### 4. `plot.Simulation()` can now warn (M9)

Given per-site intervals (from `simulate_with_shared_time_heterogeneity()` or
`join()`) and more than one site to draw, it drops the interval track and warns
`per-site intervals: pass sites=<one site> to display them`. With
`sites = <one site>` it passes that site's flat table, which is what
`RenderIntervals.js` has always expected. The js contract is unchanged.

### 5. New `PhyloPlot()` option: `download_on_click` (MIN12)

Defaults to `FALSE`. Click-to-download-the-whole-svg used to fire on any click
outside the RStudio viewer pane; it is now opt-in and gated in `PhyloPlot.js` on
`plot.options.download_on_click`. Worth a line in `man/` alongside the other
widget options.

### 6. Still open after this batch

- `toJSON(..., digits = NA)` fixes the R half of MIN8, but the C++
  newick-to-json writer (`Palantir::Phylogeny`'s json output) emits only **6
  significant digits**, so a branch length of 0.000123456789 reaches R already
  rounded to 0.000123457. Real trees are nowhere near that resolution, so this
  is cosmetic, but the payload is not exact end to end. Whoever owns the C++
  json writer could raise its precision.
- `tests_review/MIN10.R`, `MIN11.R` and `MIN12.R` assert on the INSTALLED
  javascript (browser behaviour cannot be observed from R), following the
  precedent of `MIN2`-`MIN4`. The markers they grep for are `.remove()` before
  `plot.container = d3.select` in `AlignmentPlot.js`, the absence of `n_sites`
  and `var site =` in `PhyloPlot/Tooltip.js`, and
  `plot.options.download_on_click` before `download_svg` in `PhyloPlot.js`. If
  those files are reformatted, update the markers.
- The README's `sites = 1:10` example is still wrong (sites are 0-based); not
  this agent's file.
