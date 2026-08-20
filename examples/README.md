# PalantiR command-line examples

Stand-alone R scripts that drive PalantiR from the command line. Each writes a
complete set of outputs (codon and amino-acid FASTA alignments, the full
substitution history, an interactive substitution-history plot, and a
`run_info.txt` recording parameters and versions) to an output directory, and
each accepts `--option=value` arguments; run with `--help` for the full list.
PalantiR must be installed; default trees and fitness profiles ship with the
package.

## simulate_site_heterogeneous.R

Time-homogeneous mutation-selection simulation with site-heterogeneous
fitness: each site is assigned one of a set of amino-acid fitness profiles
(by default the 50 profiles in the packaged `psi_c50_1` data set), and sites
sharing a profile are simulated as a block and combined with `join()`. Also
writes `site_profiles.tsv`, the profile assigned to each site.

```
Rscript simulate_site_heterogeneous.R --sites=200 --population-size=10000 --out=run1
```

## simulate_time_heterogeneous.R

Time-heterogeneous mutation-selection simulation: a mode tree assigns each
branch a model index and each mode runs its own effective population size
(default: the documentation example, a mammalian phylogeny with a larger
population size on the primate clade). Also writes `branch_modes.tsv`, the
mode assigned to each branch.

```
Rscript simulate_time_heterogeneous.R --population-sizes=5000,8000 --sites=300
```
