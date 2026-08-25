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

## GY94 history-level dN/dS diagnostic

`gy94_deep_tree_history_dnds.R` simulates 5,000 codon sites on a fully
bifurcating 50-taxon tree of depth 20 in synonymous-scaled branch units. Its
two 25-taxon subtrees use omega values 0.1 and 1.0. The script estimates dN/dS
from the complete substitution histories using Goldman and Yang's neutral
synonymous and nonsynonymous opportunity counts, rather than the raw event
count ratio. It writes the alignment, compressed history, tree and mode tree,
per-site counts, opportunity calculation, summary, and checksummed manifest to
a persistent artifact directory:

```sh
Rscript validation/gy94_deep_tree_history_dnds.R \
  --out=runs/gy94-deep-tree-dnds-20260825
```

The executed Jupyter notebook
[`gy94_deep_tree_history_dnds.ipynb`](gy94_deep_tree_history_dnds.ipynb)
recomputes the headline values from all 762,540 history rows, presents the
Goldman--Yang opportunity calculation, and includes three figures: assigned
versus recovered dN/dS, a complete representative-site history on the tree,
and the across-site distribution of synonymous and nonsynonymous counts. The
polished HTML reader is available both
[locally](gy94_deep_tree_history_dnds.html) and on the
[documentation site](https://dekoning-lab.github.io/PalantiR/GY94_History_Validation.html).

Rebuild the notebook source with `build_gy94_deep_tree_history_notebook.py`,
execute it top-to-bottom, and regenerate the reader with Pandoc:

```sh
python3 validation/build_gy94_deep_tree_history_notebook.py
python3 -m jupyter nbconvert --execute --to notebook --inplace \
  validation/gy94_deep_tree_history_dnds.ipynb
pandoc validation/gy94_deep_tree_history_dnds.ipynb \
  --from=ipynb --to=html5 --standalone --embed-resources --toc --mathml \
  --css=validation/notebook.css \
  --output=validation/gy94_deep_tree_history_dnds.html
```
