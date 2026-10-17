# Independent sparse-operator and preprocessing review

Review baseline: `62412f258465829d3f7110e1050d739374103781` (merged 1.0.0 candidate).
Reviewer scope: `R/preprocess.R`, input conversion in `R/utils.R`, prediction,
score diagnostics, and how the shared preprocessing object reaches a path.
The production source was not changed by this review. A focused new test file,
`tests/testthat/test-sparse-compression.R`, adds 3 blocks / 87 assertions.

## Mathematical derivation, independently of test outcomes

Write the processed data as `Z = (X - 1 mu') D^-1`, where `mu = 0` if
centering is disabled, and `D = I` if scaling is disabled. The intended
forward product is `Z B`; the transpose product is `Z' Y`.

For a fully stored column, subtracting its mean from every stored value is
an explicit dense-equivalent centering operation without adding stored
entries. For an incomplete column, its mean must also affect implicit zeros.
The implementation therefore divides the means into an explicitly applied
vector `e` (nonzero only for fully stored columns) and a remaining shift
`s = D^-1 (mu - e)`. It constructs the sparse matrix
`W = (X - 1 e') D^-1`. The two closures then evaluate:

- forward: `W B - 1 (s' B)`;
- transpose: `W' Y - s (1' Y)`.

Expanding either expression gives the intended processed-data product. The
pair are exact adjoints in real arithmetic. This also shows why adjusting
only forward products, or scaling `B` without scaling the centering shift,
would have been wrong; neither error is present.

For column j with m_j stored entries among n rows, its centered squared
norm is the sum of `(stored_x - mu_j)^2` over those entries, plus
`(n - m_j) mu_j^2` for implicit zeros. Explicit stored zeros are correctly
counted in the first term rather than the second. This formulation avoids
the cancellation-prone difference `sum(x^2) - n * mean(x)^2`.

General-real-CSC canonicalization occurs before this calculation. Expanding
symmetric/triangular/unit-diagonal representations is important because
otherwise a stored-position count would not equal a logical data-position
count. The code performs that expansion. Missing rows are found after it,
so an NA in a symmetric off-diagonal position affects both corresponding
logical observations.

The sparse randomized sketch returns `Q' Z`, expressed as the transpose of
the adjoint product. Alternating QR normalization in the power iterations
stabilizes both the observation-space and feature-space factors. The later
SVD reduction has Gram matrix `Z' Q Q' Z`, exactly the sketched model that
the solver documents. The returned scores are subsequently evaluated with
full processed Z. A fit optimizing the sketch is therefore not a certificate
that the full-data objective is globally minimized; this distinction is
already documented.

The prepared object is reused unchanged across regularization-path fits.
Only each fitted B is passed forward. There is no mutation of cached data,
scaling, means, or sketch during fitting in the reviewed code. New-data
prediction uses the training means/scales, including when the number of
new rows differs; it does not recalculate them from new observations.

## Confirmed findings

No new production correctness bug was confirmed in the sparse operator,
input canonicalization, centering/scaling, sparse sketch, or prediction
paths within this review.

A meaningful test gap was confirmed: the existing dense/sparse fitting
comparison used 12 features, k = 3, and o = 20, so the sketch had full width.
It could not establish correctness of actual sparse dimensional compression.
The sparse power-iteration implementation also lacked its own quality
regression (the existing q regression exercised the dense backend).

The new test file fills that gap with full-input reconstruction errors and
independent mathematical references, without comparing arbitrary bases or
signs. It also compares mixed large-offset/zero/tiny columns individually,
so a large column cannot mask a centering error in another column.

## Executed independent experiments

Runtime: R 4.6.1, Matrix 1.7-6, rsvd 1.0.6; BLAS/OpenMP threads fixed at one.
Raw scripts and logs are beside this report.

1. Manual explicit preprocessing references agree with both operator
   products, energy, stored means and scales for mixed sparse/constant/tiny
   columns, explicit stored zeros, triangular unit diagonals with missing
   entries, duplicated triplet positions, and fully stored large offsets.
   Input Matrix objects remain unchanged after preprocessing.
2. Dense and CSC/triplet/row-compressed zero-row prediction returns a dense
   0-by-k matrix. Duplicate feature names are rejected when matching names.
   A single NA or empty feature name is currently accepted if both sides
   contain the same name set; see limitations below.
3. A 120-by-72, exactly rank-four matrix with 1.67% nonzeros is reconstructed
   to roughly 1e-15 relative Frobenius error with k + o = 7. This holds for
   centered/uncentered and scaled/unscaled inputs and q = 0 or 2.
4. A full-rank sparse diagonal matrix with a geometric spectrum is genuinely
   compressed from 72 feature dimensions to seven sketch rows. Across five
   seeds and both centering choices, q = 4 gives squared residual at most
   1.01084 times the optimal rank-four SVD tail. q = 0 residual ratios range
   roughly 1.11 to 1.25; every probed seed benefits from refinement. The unit
   test retains three seeds and a 1.03 upper bound, with an explicit margin
   for numerical backends and no fragile factor comparisons.
5. Dense and sparse versions of that compressed fit give equal reconstruction
   quality under all centering/scaling flags. Exact reconstruction is not
   required for the full-rank case: its nonzero tail is the correct target.
6. Independent extreme-scale probes show consistent dense/sparse behavior:
   scale = TRUE handles the tested values through 1e153; a mixed-sign 1e154
   case rejects nonfinite squared scales with an explicit rescaling error.
   Subnormal data produce subnormal unscaled energy and inherit the documented
   replacement of scale factors below 1e-8 by one.

All 87 added assertions passed in targeted runs on both R 4.6.1 / Matrix
1.7-6 / rsvd 1.0.6 and R 4.0.0 / Matrix 1.3-0 / rsvd 1.0.5 / testthat
3.1.6, with no failures, errors, warnings, or skips. Existing installed
candidate libraries supplied the unchanged sparse implementation. These
are targeted test-file results, not a fresh claim that every supported
platform or a rebuilt archive was checked.

## Limits and possible future improvements

- Sparse forward/transpose products subtract a rank-one offset for incomplete
  columns. There can be ordinary floating-point cancellation, especially with
  nearly dense columns and unusual shifts. Fully stored columns are explicitly
  centered to prevent the important severe-offset case. The independent mixed
  full/incomplete test confirms that branch distinction. No failing practical
  counterexample was found here.
- Squaring before scaling can overflow even when a hypothetically standardized
  result would be finite. The implementation fails clearly and requests input
  rescaling; it does not silently claim success. A scaled sum-of-squares kernel
  could widen the representable range but would need consistent dense/sparse
  implementations and careful performance assessment. It is not required to
  fix the confirmed compressed-test coverage gap.
- The absolute small-scale cutoff 1e-8 is an established documented behavior.
  Consequently arbitrary rescaling of tiny columns is not an invariance the
  public API promises. Changing that convention should be a separate design
  decision rather than a covert numerical cleanup.
- The current name contract checks uniqueness and equal name sets. It accepts
  one missing or empty name. Requiring complete, nonempty names when aligning
  would be defensible, but would change accepted inputs and is not a confirmed
  value/dimension error under the present contract.
- Sparse input still produces dense sketch/factor/score matrices. Complexity
  is governed partly by n*(k+o) and p*(k+o), and a large requested rank removes
  the memory benefit. Deterministic and robust solvers intentionally use dense
  data; robust S is dense by its current return contract.
- Computing column means when center = FALSE is avoidable preprocessing work
  in the sparse path. It is a possible small optimization, not a measured
  performance defect. The main sparse memory-saving design is already sound.
- The QR-adjusted score convention uses a numerical rank threshold. Exact
  dependence and order-dependent allocation of energy deserve their own
  closed-form tests, owned by the independent unit-test reviewer. Neither
  adjusted energies nor raw score variance is a convergence certificate.

Recommendation: retain the current sparse algorithm and public contract,
merge the targeted compressed-quality tests, and assess any new production
changes against separately demonstrated problems rather than adding a new
backend or broad option surface in this pass.
