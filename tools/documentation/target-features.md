# Search for an approximate number of measurements

The README's feature-count example starts with two score coordinates using
three of the four iris measurements, then repeats the search for two
measurements. A measurement is counted once when any of its loading
coefficients is nonzero. The example therefore concerns the measurements
needed to calculate the complete representation.

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

## Repeat with two measurements

Petal width's small coefficient suggests investigating a smaller
representation. The second fit changes only `n_features` from 3 to 2.
It starts a separate deterministic search with the same data, preprocessing,
component count, ridge penalty, and convergence settings.

Both searches attain the requested count with converged fits. The
two-measurement search chooses sepal width and petal length after 15
candidate fits, at relative `alpha = 0.270282658390516`. Its selected fit
converges in 76 optimization iterations, with one nonzero coefficient in
each component.

| Requested and selected measurements | Relative `alpha` | Training relative reconstruction error |
| ---: | ---: | ---: |
| 3 | 0.237137370566165 | 0.561699821364799 |
| 2 | 0.270282658390516 | 0.630145039902561 |

The reconstruction error rises from 0.562 to 0.630. The selected penalty
is larger, and the retained coefficients also shrink: petal length's
magnitude changes from 0.991 to 0.901, while sepal width's changes from
0.233 to 0.104. Raw component signs differ between fits, so magnitudes or
the saved sign-oriented coefficients are the useful comparison.

These are complete refits under different selected penalties. The first
13 penalties and their initial warm starts are the same, but the requested
count changes which candidate initializes the refinements. The complete
search records preserve those differences. The errors compare the returned
penalized models; they do not measure the isolated contribution of the
three-measurement fit's small petal-width term or the best reconstruction
obtainable with two measurements. A choice for future observations needs
held-out assessment. Candidate and iteration counts are recorded for
reproducibility and are not timing measurements.

## Reproduce the example

From the repository root, with the package dependencies and development
tool `pkgload` available:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript --vanilla tools/documentation/target-features.R
```

The script checks the selected loading count and score calculation for
both fits, saves both searches, and draws the three-measurement PNG and
SVG. The result snapshot retains the input, parameters, chosen fits,
search histories, and source-file
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

The snapshot has format version 2. Its `fits$three` and `fits$two` entries
contain the respective parameters, target results, selected variable names,
display coefficients, and reconstruction errors. Each CSV identifies its
search with `requested_features`. Evaluation and warm-start indices belong
to that search, so they restart when the requested count changes.

## Saved evidence

| File | Contents |
| --- | --- |
| [target-features-iris.rds](results/target-features-iris.rds) | Input matrix, both sets of parameters, complete selected fits, search records, and generation metadata. |
| [target-features-summary.csv](results/target-features-summary.csv) | One row per request with achieved count, selected variables and penalty, iteration counts, and relative reconstruction error. |
| [target-features-search.csv](results/target-features-search.csv) | Every penalty, signed count difference, convergence result, and source of the warm start, grouped by requested count. |
| [target-features-loadings.csv](results/target-features-loadings.csv) | Raw and displayed coefficients, labeled by requested count, measurement, and component. |
| [target-features-session-info.txt](results/target-features-session-info.txt) | Environment, thread settings, and source identity. |

The figure is stored as
[PNG](../../man/figures/target-features-iris.png) and
[SVG](../../man/figures/target-features-iris.svg).
The figure retains the three-measurement example. Its left panel follows
the actual evaluation order, while the right panel shows which measurements
enter each score. A horizontal reference marks
the requested count. No interpolation of the penalty-to-support relation
is assumed by the connecting line, which only links consecutive trials.

## Interpretation

The example demonstrates two achieved counts for one dataset under the
stated settings. The search selects from converged candidates by
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
