# Engine validation notebook

`engine_validation.Rmd` tests the core simulation engine against exact
expectations and renders to a standalone HTML report with a PASS/FAIL
verdict per test:

1. **Branch-length calibration** — one unit of branch length delivers one
   expected substitution of the scaled class (synonymous, non-synonymous, or
   any) per codon site.
2. **Stationary distribution** — tip state frequencies match the model's
   computed equilibrium (chi-square and sequence logos).
3. **Transient between two equilibria** — in a time-heterogeneous simulation
   the state distribution along the branch follows the master equation of the
   new model, integrated in the branch's own units, under both rescale
   methods. This test found the defect fixed in version 1.2.1.

Rerun with a different seed or size:

```r
rmarkdown::render("engine_validation.Rmd", params = list(seed = 1, n_sites = 3000))
```

Requires PalantiR and rmarkdown; `ggseqlogo` is optional (barplots are drawn
without it).
