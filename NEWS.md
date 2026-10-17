# sparsepca 1.0.0

## Correctness and observable behavior

- Correct `robspca()`'s objective to use the complete residual
  `X - X %*% B %*% t(A) - S`. Previously the sparse error matrix was
  charged a penalty while omitted from the reconstruction loss, which
  could increase the recorded objective and terminate a fit prematurely.
  Thanks to Maria Sans-Fuentes for the report and correction. Robust fits
  with nonzero sparse residuals can now run longer and produce different
  fitted values; the intended optimization model is unchanged.
- Honor the `q` argument to `rspca()`; previous versions always used two
  power iterations. Calls with the default `q = 2` retain that setting.
- Make rank-one, zero-data, rank-deficient, and single-variable cases
  well-defined, retain matrix dimensions, and validate unsupported inputs
  and invalid scalar arguments before numerical operations.
- Correct `eigenvalues` and `sdev` to report sequentially adjusted score
  energy and its square root. The previous values came from an auxiliary
  Procrustes SVD and were not the variances of the returned sparse scores.
  These field names and the summary matrix layout are retained; their
  numerical values can change for penalized fits. `score_variance` records
  marginal score energy separately. With centered data these quantities
  are adjusted and marginal sample variances, respectively.
- Compute `var` as the total squared energy of the processed data divided
  by `n - 1`. This is total sample variance when centering is enabled.
  For `center = FALSE`, it now consistently includes the column means,
  matching the uncentered score-energy convention.
- Retain all existing function names, argument positions and defaults,
  result classes, and existing result components. Sparse loading vectors
  remain in their fitted scale; they are not normalized after optimization.

## New capabilities

- Add `penalty = "l0"` for hard thresholding, with the existing ridge
  penalty providing the ZeroNet variant when `beta > 0`. The default
  remains elastic net. The new arguments follow all existing arguments.
- Add `B_init` to initialize a fit from a finite loading matrix, and
  `spca_path()` to fit a sequence of penalties using a shared matrix
  factorization or randomized sketch. Paths retain the supplied order and
  report per-component nonzero counts and the union of selected variables.
  Penalty paths do not impose an exact requested cardinality.
- Add `predict()` methods for all three fitted classes. Predictions use
  stored centering and scaling, accept numeric data frames and Matrix
  inputs, and align a complete set of named variables with the training
  columns. Empty dimension names are normalized across supported Matrix
  versions. Calling `predict()` without new data returns the training scores.
- Support sparse Matrix input in `rspca()` with implicit centering and
  scaling. Dense sketches, loadings, and scores are still needed;
  `spca()` and `robspca()` explicitly use dense working matrices.
- Return iteration, objective-convergence, and effective-penalty
  diagnostics so users can inspect a fit instead of inferring convergence
  from printed output.

## Maintenance and documentation

- Extend the independent numerical tests to cover genuinely compressed
  sparse sketches, finite-difference reduced-objective gradients, correlated
  loading optimality conditions, active robust Huber optima, changes of
  measurement units, and numerical safeguard branches. Fault-injection
  experiments verify that the new safeguards detect their intended failures.
- Clarify that ordinary-PCA recovery at zero penalties assumes the default
  SVD initialization, and distinguish objective convergence from coefficient
  accuracy or first-order optimality of every returned factor.

- Avoid computing unused left singular vectors during deterministic
  initialization. Reuse thin matrix products during optimization instead
  of constructing a full reconstruction residual on every iteration.
  Near an exact fit, evaluate the residual directly in bounded blocks
  to avoid cancellation in differences of large squared norms.
- Preallocate objective histories, share common solver and preprocessing
  code, and remove duplicate products from robust updates. These changes
  preserve the intended optimization updates while reducing allocation
  and repeated work.
- Expand every help topic and add documentation for prediction,
  regularization paths, parameter scaling, reconstruction, and the
  interpretation of correlated sparse component scores.
- Expand the README into a self-contained tutorial with the two-factor
  model, parameter guidance, held-out prediction, and reproducible iris,
  grouped-measurement, and corrupted-sensor applications. Restore all five
  historical eigenface panels with provenance and interpretation, and
  retain the measured performance tables and figures.
- Deepen all seven help topics with practical guidance and independent
  examples. Check labeled README chunks and help examples against the
  installed public API, fail unexpected warnings, and record explicit
  installation skips and source provenance. Save new application results
  separately from their figure rendering.
- Add a package citation for the final 2020 SIAM publication and update
  author and maintainer metadata.
- Add numerical regression tests, package-check workflows, and
  reproducible benchmark and release-validation scripts. Development
  measurements and README assets are excluded from the source archive.

## Dependency policy

Runtime code uses Matrix's established sparse-matrix interfaces, base R
linear algebra, and rsvd's randomized range finder. Minimum versions are
R 4.0.0, Matrix 1.3-0, and rsvd 1.0.5. The R floor accommodates the
dependency toolchain without requiring the newest R release. Matrix does
not need to be attached by the user. The test suite requires testthat
3.0.0 or later.
