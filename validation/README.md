# Engine validation notebook

`engine_validation.Rmd` tests the core simulation engine against exact
expectations and renders to a standalone HTML report with a PASS/FAIL
verdict per test:

1. **Event-count branch-length calibration** — one unit delivers one expected
   synonymous, nonsynonymous, or total substitution per codon under the
   corresponding event-count gauge.
2. **dS branch-length calibration** — synonymous histories recover input dS
   after division by synonymous opportunity, for homogeneous and
   time-heterogeneous GY94 and mutation--selection models.
3. **Stationary distribution** — tip state frequencies match the model's
   computed equilibrium (Holm-adjusted chi-square comparisons and sequence
   logos).
4. **Transient between two equilibria** — in a time-heterogeneous simulation
   the state distribution along the branch follows the master equation of the
   new model, integrated in the branch's own units, under both rescale
   methods; its checkpoint comparisons use Holm adjustment. This test found
   the defect fixed in version 1.2.1.

Rerun with a different seed or size:

```r
rmarkdown::render("engine_validation.Rmd", params = list(seed = 1, n_sites = 3000))
```

Requires PalantiR and rmarkdown; `ggseqlogo` is optional (barplots are drawn
without it).

## GY94 history-level dN/dS diagnostic

`gy94_deep_tree_history_dnds.R` simulates 5,000 codon sites on a fully
bifurcating 50-taxon tree of depth 20 using the historical
`"synonymous-per-codon"` event-count gauge. Its two 25-taxon subtrees use
omega values 0.1 and 1.0. The script estimates dN/dS
from the complete substitution histories using Goldman and Yang's neutral
synonymous and nonsynonymous opportunity counts, rather than the raw event
count ratio. It writes the alignment, compressed history, tree and mode tree,
per-site counts, opportunity calculation, summary, and checksummed manifest to
a persistent artifact directory:

```sh
Rscript validation/gy94_deep_tree_history_dnds.R \
  --out=runs/gy94-deep-tree-dnds-20260825
```

The primary feature demonstration is the executed R-kernel notebook
[`gy94_deep_tree_history_dnds_R.ipynb`](gy94_deep_tree_history_dnds_R.ipynb).
It loads the PalantiR checkout, constructs the tree and branch-heterogeneous
GY94 model, performs the complete 5,000-site simulation, examines the resulting
R objects through PalantiR's native interactive alignment and history widgets,
and reconstructs the headline values from all 762,540 events. The full
fixed-seed simulation takes only a few seconds, so the notebook does not use a
reduced toy run.

Install and register an R kernel once if Jupyter does not already list one:

```r
install.packages(c("IRkernel", "pkgload"))
IRkernel::installspec()
```

Rebuild and execute the live notebook from the repository root:

```sh
python3 validation/build_gy94_deep_tree_history_r_notebook.py
python3 -m jupyter nbconvert --execute --to notebook --inplace \
  --ExecutePreprocessor.kernel_name=ir \
  --ExecutePreprocessor.timeout=240 \
  validation/gy94_deep_tree_history_dnds_R.ipynb
```

The separately preserved Python notebook
[`gy94_deep_tree_history_dnds.ipynb`](gy94_deep_tree_history_dnds.ipynb)
is an independent audit: it streams the committed production history rather
than invoking the simulator and checks the same event totals and Goldman--Yang
normalization. It also provides three static figures: assigned versus recovered
dN/dS, a complete representative-site history, and the across-site distribution
of synonymous and nonsynonymous counts. Its polished HTML reader is available
both
[locally](gy94_deep_tree_history_dnds.html) and on the
[documentation site](https://dekoning-lab.github.io/PalantiR/GY94_History_Validation.html).

Rebuild the independent audit and regenerate its reader with Pandoc:

```sh
python3 validation/build_gy94_deep_tree_history_notebook.py
python3 -m jupyter nbconvert --execute --to notebook --inplace \
  validation/gy94_deep_tree_history_dnds.ipynb
pandoc validation/gy94_deep_tree_history_dnds.ipynb \
  --from=ipynb --to=html5 --standalone --embed-resources --toc --mathml \
  --css=validation/notebook.css \
  --output=validation/gy94_deep_tree_history_dnds.html
```
