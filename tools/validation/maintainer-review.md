# Mathematical and public-contract review for sparsepca 1.0.0

Review date: 2026-10-17. This document records implementation decisions;
the release evidence report records the final source/archive identities and
which checks actually ran.

## Starting point and scope

The GitHub baseline is `fa896289a68801a4798f53948ff6c2b8f7c54e1c`
(DESCRIPTION version 0.1.0). The R files in that revision, the published
CRAN 0.1.2 archive, and the maintainer-supplied 0.1.2 ZIP were compared
byte for byte and are identical. The published archive is therefore the
appropriate downstream baseline. Metadata differs: 0.1.2 adds the arXiv
reference and CRAN publication/build fields.

The package has three S3 generics, their default methods, and print/summary
methods: `spca`, `rspca`, and `robspca`. There are no compiled routines or
bundled datasets. Help pages are generated from roxygen source. Previously
the three test files were never invoked by `R CMD check`: there was no
`tests/testthat.R`, and the files attempted `library(spca)` instead of the
actual package name `sparsepca`.

The intended model is variable-projection sparse PCA. It distinguishes
the sparse coefficient matrix B from an orthonormal reconstruction factor
A. Loadings are penalized coefficients, with no normalization after
optimization. Scores are X B for the preprocessed input. They can be
correlated, so ordinary PCA identities cannot be transferred without
checking their assumptions.

## Prioritized findings

| Priority | Finding | Evidence and resolution |
| --- | --- | --- |
| Confirmed correctness bug | Robust objective omits S from squared loss | Paper equation 17 and direct objective evaluation agree on `0.5 * sum((R-S)^2)`. A controlled contaminated example stops after two old iterations; correcting only this term continues toward a lower objective. Regression independently reconstructs the objective from returned factors. |
| Confirmed correctness bug | `rspca(q=...)` ignored | The old call hardcodes `q=2`. Fixed-seed old fits at q=0 and q=4 are identical. The new quality test uses a slowly decaying spectrum and verifies the effect of extra power iterations. |
| Confirmed reporting bug | `eigenvalues` is an intermediate Procrustes singular-value vector | With one component and ridge shrinkage, this quantity differs analytically from the variance of the returned score. Replace it with sequential QR-adjusted score energy, keeping the historical field names. |
| Confirmed reporting inconsistency | Uncentered numerator divided by centered total variance | Use processed input squared energy divided by n-1 consistently. Centered fitting still reports sample variance; uncentered fitting reports second moments with denominator n-1, matching the PCA reference convention. |
| Confirmed boundary bug | Zero data causes infinite step sizes and undefined stopping comparisons | Return finite zero scores/objective and dimensionally valid factors, with explicit convergence. Reject invalid scalar controls and unsupported nonfinite/complex input before optimization. |
| Test weakness | Signed sums can cancel reconstruction errors | Replace with Frobenius residuals, optimal spectral tails, subspace projectors, independent scalar optima, and explicitly reconstructed objectives. |
| Performance opportunity | Repeated full residuals and unused SVD vectors | Keep a thin/right SVD factor, request no unused left vectors, use products with thin loading factors, and preallocate objective history. |
| Memory opportunity | Every sparse input is densified | Randomized fitting now uses explicit sparse products with implicit centering/scaling. Deterministic and robust fitting retain their documented dense workspace. |
| Documentation gap | Coefficients, correlations, variances, reconstruction and prediction are conflated | Every help topic and the README now distinguish these quantities and show stored preprocessing, practical tuning, and memory limits. |
| Requested enhancement | Hard thresholding/ZeroNet (#2) | Add an optional cardinality penalty with a tested hard-threshold proximal step. Elastic net remains the default. |
| Requested enhancement | Warm starts and sparsity exploration (#3) | Add `B_init` and a path that reuses preprocessing/factorization/sketch. Record both per-component nonzeros and the union of selected variables. Exact target-cardinality optimization is not claimed. |

Issue #1 is not an outstanding indexing bug: its author retracted it after
observing that the objective vector starts empty. The current indexing
compares the two latest completed iterations. Issues #4–#7 are historical
installation/interpretation/prediction discussions; documentation and new
prediction methods incorporate their useful lessons.

## Numerical implementation decisions

For deterministic fitting, let C = D V' be the right SVD factor of the
processed input. It has the same Gram matrix as X. For randomized fitting,
the same construction starts from the small sketch Q' X. The orthogonal
update is `polar(crossprod(C, C %*% B))` and the negative smooth gradient
is `crossprod(C, C %*% (A-B)) - beta_effective * B`.

The existing step size and relative penalty scaling are preserved:
`alpha_effective = alpha * d1^2`, `beta_effective = beta * d1^2`, and
step size `1 / (d1^2 + beta_effective)`. Algebraically the step is computed
as `(B + crossprod(C, C %*% (A-B)) / d1^2) / (1+beta)`, avoiding an
extra large product. Soft thresholding uses `alpha/(1+beta)`; hard
thresholding uses its square root after multiplication by two. The ridge
term stays in the smooth gradient for both penalties. Independent tests
compare scalar orthogonal examples with their analytic solutions.

The squared residual is evaluated through
`sum(C^2) - 2*sum((C A)*(C B)) + sum((C B)^2)`, which is valid because
A' A = I. Near exact reconstruction, subtraction loses relative accuracy.
When the result is at most `1e-8 * sum(C^2)`, bounded row blocks instead
evaluate the actual residual C - (C B) A'. This is a numerical fallback,
not clipping an inaccurate loss to zero. The objective-increase check
allows floating-point uncertainty relative to both the objective and
data energy, while significant increases warn and do not report
convergence. The tolerance controls objective stagnation, not a proof of
stationarity or global optimality.

Robust fitting uses `crossprod(X-S, X B)` for the orthogonal update and
`crossprod(X, (X-S) A - X B)` for the gradient. This avoids constructing a
second full gradient residual. The sparse-error update remains soft
thresholding of X - (X B) A'. Its objective now charges the reconstruction
loss after subtracting S, plus the loading and sparse-error penalties.
S remains a dense n-by-p R matrix, as in the public contract.

## Sparse products and numerical stability

Structural, pattern, logical, triplet and row-compressed sparse Matrix
objects are converted to a real general column-compressed representation
through documented coercion classes. Symmetric entries and unit diagonals
are expanded correctly in sparse storage. The implementation uses explicit
`Matrix::crossprod` dispatch, avoiding version-sensitive assumptions about
base product generics.

Forward and transpose products implement `(X-1*mu') D^-1` without storing
the full centered input. Column energy is accumulated from deviations of
stored values plus contributions from implicit zeros; it does not subtract
`n*mean^2` from `sum(x^2)`. Fully stored columns are centered directly in
their sparse values to avoid cancellation from a large constant offset.
Their means follow base `colMeans` exactly. Tests first check both products
and adjoint consistency against explicit dense matrices, then compare
fits, and finally profile allocations to reject an n-by-p dense array.

The method is sparse in input storage, not in every intermediate. The
basis, sketch, coefficient matrices and scores are dense. Matrix size,
sketch width and requested rank still govern memory requirements.

## Compatibility decisions

- All existing argument positions/defaults, generics, classes and result
  components are retained. `penalty` and `B_init` are appended. The original
  print return value (invisible rounded loading matrix) and four-row summary
  layout are preserved.
- `loadings` is not renormalized. Sign and basis ambiguity are handled in
  tests using reconstruction/subspaces, not byte-identical factors.
- Corrected variance fields are deliberately observable. QR adjustment is
  sequential in returned component order, with numerically dependent
  directions assigned zero. It is not a covariance eigen-decomposition.
  Marginal score energy is available separately as `score_variance`.
- Missing training rows are still omitted with a warning. At least two
  complete rows are required. Constant scale factors below `1e-8` still
  become one. Rank requests above min(n,p) are still clipped.
- Randomized fits use the caller's R RNG. No public method resets the
  seed. A randomized path draws one sketch and reuses it across penalties.
- Prediction is additive API: it uses stored means/scales and original
  coefficient ordering, including a validated reordering of unique named
  columns. It does not refit outliers or silently remove new observations.
- Runtime minima are R 4.0.0, Matrix 1.3-0 and rsvd 1.0.5, with exact
  lower-bound testing separate from current dependency versions.

## References reviewed

- Erichson et al. (2020), *Sparse Principal Component Analysis via Variable
  Projection*, SIAM Journal on Applied Mathematics 80(2), 977–1002.
  <https://doi.org/10.1137/18M1211350>;
  <https://arxiv.org/html/1804.00341v3> (robust objective and Appendix A).
- Zou, Hastie and Tibshirani (2006), *Sparse Principal Component Analysis*,
  Journal of Computational and Graphical Statistics 15(2), 265–286.
  <https://doi.org/10.1198/106186006X113430> (adjusted score variance).
- R Core, *Writing R Extensions* and CRAN Repository Policy:
  <https://cran.r-project.org/doc/manuals/r-release/R-exts.html>;
  <https://cran.r-project.org/web/packages/policies.html>.
- The revised rSVD repository was reviewed for preprocessing, documentation,
  and source-archive validation conventions: <https://github.com/erichson/rSVD>.

No claim of universal robustness follows from this review. Full external
validation, baseline/candidate comparisons, benchmark measurements and
remaining environment limitations belong in the accompanying release report.
