# Search for an approximate number of measurements

The README's feature-count example asks for two score coordinates using
three of the four iris measurements. A measurement is counted once when
any of its loading coefficients is nonzero. The example therefore concerns
the measurements needed to calculate the complete representation.

The input is the four numeric columns of R's
[iris dataset](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/iris.html),
containing 150 observations. Species labels are unused. The fit uses
`spca_target()` with `k = 2`, `n_features = 3`, `beta = 0.001`, centering
and scaling, `tol = 1e-8`, and `max_iter = 3000`. The default search starts
with at most 13 relative penalties, allows 30 total fits, and stops at the
first converged candidate with the requested count.

The initial grid and any refinements share one data transformation and one
factorization. Each initial candidate starts from the preceding loadings;
refinement starts from the closest converged candidate so far. Intervals
whose observed endpoint counts fall on either side of the request are
investigated first. The search can examine several such intervals and does
not assume that support changes monotonically with the penalty.

The saved run found three selected measurements after 16 candidate fits,
at relative `alpha = 0.237137370566165`. The selected fit converged in
195 optimization iterations. Its first component uses petal length and
petal width; the second uses sepal width. One earlier candidate reached
the 3,000-iteration limit and was excluded from selection. The search
record retains that candidate and its warm-start relationship.

The relative reconstruction error on the standardized input is
0.561699821364799, calculated as the Frobenius norm of the reconstruction
residual divided by the Frobenius norm of that input. This quantity is a
training-data fit diagnostic. Petal width's coefficient is approximately
0.000503, compared with 0.991 for petal length in the same component.
The requested count is attained, while the coefficient magnitudes show
that those two measurements contribute very differently.

## Reproduce the example

From the repository root, with the package dependencies and development
tool `pkgload` available:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript --vanilla tools/documentation/target-features.R
```

The script checks the selected loading count and score calculation, saves
the search result, and draws the PNG and SVG. The result snapshot retains
the input, parameters, chosen fit, search history, and source-file
checksums. Source commit and working-tree state are recorded together so
an uncommitted candidate is not mistaken for the committed source alone.
The session record identifies R, dependencies, BLAS/LAPACK, and thread
settings.

To change the presentation without fitting again:

```sh
Rscript --vanilla tools/documentation/target-features.R --plot-only
```

This reads the saved result and uses base R graphics. Loading columns are
sign-oriented for display so their largest-magnitude coefficient is
positive; both raw and displayed coefficients are retained. This display
choice changes neither selected variables nor the saved fit.

## Saved evidence

| File | Contents |
| --- | --- |
| [target-features-iris.rds](results/target-features-iris.rds) | Input matrix, parameters, complete selected fit, search record, and generation metadata. |
| [target-features-summary.csv](results/target-features-summary.csv) | Requested and achieved counts, selected penalty, iteration counts, and relative reconstruction error. |
| [target-features-search.csv](results/target-features-search.csv) | Every penalty, signed count difference, convergence result, and source of the warm start. |
| [target-features-loadings.csv](results/target-features-loadings.csv) | Raw and displayed coefficients, labeled by measurement and component. |
| [target-features-session-info.txt](results/target-features-session-info.txt) | Environment, thread settings, and source identity. |

The figure is stored as
[PNG](../../man/figures/target-features-iris.png) and
[SVG](../../man/figures/target-features-iris.svg).
Its left panel follows the actual evaluation order, while the right panel
shows which measurements enter each score. A horizontal reference marks
the requested count. No interpolation of the penalty-to-support relation
is assumed by the connecting line, which only links consecutive trials.

## Interpretation

The example demonstrates an achieved count for one dataset and one
optimization path. The search selects from converged candidates by
absolute count difference, then fewer variables, then earlier evaluation.
It does not compare penalized objective values across different penalties.
The returned coefficients solve the original penalized model; they have
not been clipped to create the requested count.

An approximate target expresses a practical preference. Other datasets or
warm starts can skip a requested count or reach different local solutions.
The result reports the achieved count and `target_met`; a search that
misses returns its closest converged candidate. Reaching the count does
not establish a unique subset, statistical significance, or the best
predictive representation. The figure and its reconstruction error are
illustrative fit diagnostics, rather than a feature-selection benchmark.
