# Second review: test adequacy and implementation soundness

Review date: 17 October 2026. Baseline: merged `sparsepca` 1.0.0 candidate
`62412f258465829d3f7110e1050d739374103781`.

## Conclusions and changes

The second review found meaningful test gaps and an overly broad documentation
claim, but no additional confirmed defect in the fitting updates. It independently
derived the numerical operations, reviewed all existing tests and public methods,
and used mathematical references rather than treating a green test suite as a
proof. The executable fitting, preprocessing, prediction, and path code is
unchanged. DESCRIPTION and NAMESPACE are unchanged as well; exact hashes are in
[source-equivalence.json](source-equivalence.json).

The suite increases from 29 blocks / 621 assertions to **40 blocks / 784
assertions**. All 784 passed under R 4.6.1 against the unchanged installed rc2
implementation. Full checks of the subsequently rebuilt source archive and
other environments are recorded in the release handoff, not inferred from this
initial source-level run. A separate run exercised all 784 assertions with
R 4.0.0, Matrix 1.3-0, rsvd 1.0.5, and the exact declared testthat 3.0.0 floor.
Compiling that old test tool on modern glibc required its built-in
`CATCH_CONFIG_NO_POSIX_SIGNALS` option; its archive and R implementation were
unchanged. The initial compile failure and an edition-selection harness pilot
are retained separately from the final passing edition-3 run.

The documentation now specifies that unpenalized ordinary-PCA recovery assumes
the default SVD initialization. It also explains the ordering of the final
A/B/S updates and distinguishes objective stagnation from coefficient error or
first-order optimality.

## Why the additional tests are useful

| Gap in the earlier suite | Added reference or property | Assertions |
|---|---|---:|
| Sparse/dense fits used p=12, k=3, o=20, hence a full-width sketch | Actual rank-four sparse input compressed to seven sketch rows; full-rank sparse spectrum tested against the optimal SVD tail and q=0 versus q=4 | 42 |
| Sparse and dense helpers shared some preprocessing code | Independent explicit centered/scaled targets with mixed large offsets, implicit zeros, tiny columns, and immutable input | 45 |
| Penalty references were primarily orthogonal scalar problems | Finite differences of the nuclear-norm reduced objective, followed by the independently computed proximal step | 2 |
| Small objective change alone did not establish correlated loading optimality | Full-residual elastic-net KKT conditions, orthogonality, and the exact Procrustes trace maximum | 6 |
| Active robust regression evaluated the objective but did not independently locate an optimum | Profiled one-dimensional Huber root and active/zero l0 comparisons; gamma=0 exact limit | 23 |
| Near-exact residual cancellation and stopping classifications were not forced | Direct residual at loss about 4.25e-17; descent, roundoff, significant increase, and zero objective | 17 |
| QR result checks partly shared the QR family of algorithms | Closed-form score geometry with zero/dependent columns followed by independent ones | 2 |
| Relative regularization and robust gamma units were not jointly challenged | Global change-of-units relation for all three methods and both penalties, including a negative multiplier | 26 |

These 163 additional assertions are focused on different failure modes. The
original tests remain useful: matched-preprocessing PCA, controlled spectral
tails, sign-invariant subspaces, true robust objectives, analytic shrinkage,
input validation, paths, prediction, RNG behavior, and sparse allocation checks.
The [complete test audit](tests/audit.md) maps all 29 earlier blocks to the
specific evidence each supplies and its limitations.

### Fault-injection checks

Seven deliberate, process-local faults were detected. No repository or installed
implementation files were modified by these experiments.

| Injected fault | Evidence that the new tests detect it |
|---|---|
| Remove the stable direct-residual fallback | The near-exact relative loss check fails |
| Treat a meaningful objective increase as convergence | The convergence-classification check fails |
| Ignore data energy in the roundoff allowance | The tiny-increase classification checks fail |
| Ignore the QR pivot when assigning adjusted energies | The known score-geometry check fails |
| Omit ridge shrinkage | Both finite-difference and correlated-KKT blocks fail |
| Reverse the smooth gradient sign | Both finite-difference and correlated-KKT blocks fail |
| Treat relative alpha as an absolute threshold | Both finite-difference and correlated-KKT blocks fail |

The mutation runs used an earlier focused geometry snapshot with seven
assertions and independently calculated KKT penalty coefficients. A parsed-
expression comparison confirms that the final version adds one orthogonality
assertion and retains all prior conditions unchanged. The saved mutation evidence is carried forward with that explicit
qualification. The full final test suite is executed independently. See
[geometry review](tests/geometry-review.md),
[mutation results](tests/mutation-summary.csv), and
[geometry mutation results](tests/geometry-mutation-summary.csv).

## Mathematical reasoning behind the assessment

### Reduced factor and orthogonal update

For the processed input, let X = U D V' and C = D V'. Then C'C = X'X.
The deterministic reduction preserves the entire penalized objective for any
A and B, not only the PCA initializer. For orthonormal A,

```math
\frac12\|C-CBA^T\|_F^2
=\frac12\|C\|_F^2+\frac12\|CB\|_F^2-\mathrm{tr}(A^TC^TCB).
```

Consequently, the polar factor of C'CB has the correct orientation and globally
solves the fixed-B orthogonal subproblem. The fast residual identity follows
from the same expansion. Evaluating the actual residual in blocks near an
exact reconstruction addresses cancellation in the difference of energies.

For randomized fitting the same statements apply to the sketch Q'X. They do
not imply that its minimizer is the exact full-data minimizer. Full-data scores
are computed afterward with the original processed-data operator.

### Loading step and its descent argument

For fixed A, the smooth loading gradient is

```math
g=C^TC(B-A)+\lambda_2B,
\qquad L=\|C\|_2^2+\lambda_2.
```

The code implements the proximal step with step length 1/L, relative penalties
lambda1=alpha*||C||2^2 and lambda2=beta*||C||2^2, and the correct soft or hard
threshold. Hard thresholding compares the scalar minimum at zero with the
active minimum, with threshold sqrt(2*lambda1/L).

The fixed-A quadratic Taylor majorizer with curvature L dominates the smooth
loss. The exact scalar proximal minimizers do not increase that majorizer,
and the exact orthogonal update also cannot increase the objective. This
establishes nonincrease in exact arithmetic for either penalty. It does not
establish a global solution or convergence of every iterate sequence.

The finite-difference test uses the independently eliminated value

```math
v(B)=\frac12\|X\|_F^2+\frac12\|XB\|_F^2-\|X^TXB\|_*.
```

Its test point is safely full column rank: the smaller singular value of
X'XB is 21.14. The observed loading-step discrepancy is about 1.87e-11,
versus a 1e-7 test allowance, with comfortably nonzero distances to threshold
ties. The correlated example also has both active and zero coefficients; its
normalized KKT residual is about 9.3e-8 at objective tolerance 1e-13.

The value function is generally nonsmooth at rank loss; at B=0 with nonzero
input this is unambiguous. The fixed-A majorizer argument remains valid there.
The implementation should not be justified by assuming that the projected
value has a globally Lipschitz gradient. The primary paper's Appendix C.3
invokes convexity of its penalty, so that particular proof is not a blanket
claim for every l0/rank-deficient case. See the
[full independent solver derivation](solver/second-pass-solver-review.md).

### Robust fitting

At fixed S, the correct cross-product for the orthogonal step is (X-S)'XB.
The loading gradient is X'(XB-(X-S)A)+lambda2*B. Soft thresholding the low-rank
residual X-XBA' gives the exact S update, and the loss uses X-XBA'-S.
These are the operations in the implementation.

Eliminating S gives the entrywise Huber loss. The new scalar oracle minimizes
that independently: an active-outlier example gives coefficient 0.5295209,
compared with fitted 0.5295218 at objective tolerance 1e-13. The objective
agrees more closely than the coefficient, as expected from objective-based
stopping. With gamma=0 and positive loading penalties, B=0,S=X is the exact
zero-objective solution; both penalty branches reach that reference case.

Robustness is specific to this model. Scores still use XB, preprocessing uses
ordinary means/scales, and sparse errors are entrywise. This does not imply
robustness to arbitrary leverage, rowwise contamination, or a guarantee of
recovering the scientifically correct support.

### Sparse operators and public methods

For Z=(X-1*mu')D^-1, the forward and transpose operators implement ZB and Z'Y
through sparse products and rank-one mean corrections. Means are directly
subtracted from fully stored columns to protect large-offset calculations.
Stored deviations and implicit-zero contributions determine squared column
norms without subtracting n*mean^2 from sum(x^2).

Canonical general CSC conversion makes structural symmetric/triangular entries
and unit diagonals explicit in sparse storage. This happens before counting
entries or omitting missing rows. Independent probes included a true symmetric
sparse off-diagonal missing value, triplet duplicates, stored zeros, and empty
prediction batches. No new value, dimension, or dispatch defect was found.
The [operator review](operators/REVIEW.md) contains the complete derivation.

Prediction uses stored preprocessing and aligns uniquely named features.
Paths reuse preprocessing/factorization or a single random sketch; they change
only the next loading initializer. These preserve the stated contracts.

## Findings that affect interpretation

1. **Default initialization matters.** For X formed from positive and negative
   diagonal entries (4,3,2), a one-component zero-penalty fit started in the
   weakest coordinate stays there: squared residual 50 rather than the optimal
   26. This is a stationary subspace, not a wrong gradient. The default SVD
   start recovers the leading PCA subspace. README and help now distinguish them.
2. **Objective tolerance is not a stationarity tolerance.** A is updated before
   the final B/S blocks. Across ten controlled contaminated problems, redoing
   only the last A update lowered the objective by at most 1.69e-6 relative for
   ordinary fits and 2.16e-6 for robust fits at tol=1e-5. Normalized A residuals
   reached about 0.0012. At tol=1e-10, objective gaps fell below 4e-11 and those
   residuals below 5.4e-6. This is a documented accuracy distinction, not a reason
   to silently change iteration order or add a polishing cycle.
3. **Finite arithmetic has limits.** Squared norms can overflow at extreme
   values even when a hypothetical standardized result is finite. The package
   requests rescaling rather than silently returning invalid results. The
   established absolute scale cutoff 1e-8 is deliberately retained.
4. **The randomized optimum and full-data optimum differ.** Convergence of
   optimization on the sketch does not establish approximation accuracy.
   Oversampling, power iterations, and full-input quality assessment remain
   separate decisions.

## Additional functionality and optimizations

No new public fitting option or production optimization was justified by a
confirmed defect in this pass. The first update already removed the major
redundant residuals/unused SVD vectors and introduced factorization reuse and
sparse input. Keeping those reviewed paths stable is preferable to adding an
unmeasured backend or more tuning controls.

Prioritized follow-on possibilities:

| Proposal | Benefit | Requirement before implementation |
|---|---|---|
| Optional first-order/proximal diagnostics | Distinguish objective stagnation from block optimality | Define A, B, and S residuals, l0 fixed-point interpretation, and sketch versus full-data meaning; do not change default stopping silently |
| Reconstruction convenience method | Centralize inverse preprocessing and optional robust error addition | Define original-unit/processed-unit behavior and dense allocation explicitly; current documented formula already works |
| Lean regularization-path results | Reduce repeated n-by-k score storage for many penalties | Demonstrate a large-path workload and preserve default return/prediction behavior |
| Small-p cached Gram calculation | Potentially reduce repeated thin products | Measure an actual speed benefit and the extra p-by-p memory; avoid compromising stability |
| Skip unused sparse means when center=FALSE | Removes a preprocessing pass | Establish meaningful end-to-end benefit; currently an unmeasured micro-optimization |
| Fully sparse/blocked robust workspace | Reduce memory for large robust problems | Separate algorithm/workspace project with explicit dense-output cost and independent equivalence tests |

Multistart, exact-cardinality constraints, adaptive sketch sizing, and an
accelerated optimizer are larger additions. None is necessary to make the
present routines mathematically consistent; each needs its own objective,
behavioral contract, reference tests, and measurements. Exact cardinality
remains the unresolved part of issue #3.

## Reproducing the evidence

Run the package's normal tests after installing it. The manual fault-injection
scripts are excluded from the CRAN archive and should run in fresh R processes:

```sh
Rscript tools/validation/second-review/tests/mutation-checks.R . /tmp/sparsepca-safeguard-mutations
Rscript tools/validation/second-review/tests/geometry-mutation-checks.R . /tmp/sparsepca-geometry-mutations
```

Supply R_LIBS_USER for the desired installed package/dependency combination.
Scripts accept output directories and do not alter installed files. The saved
raw CSVs describe the original executed experiments; rerunning the geometry
script now exercises the strengthened final block. Versioned text logs have trailing whitespace normalized. Original harnesses,
unmodified logs, and full RDS objects are preserved separately in the review evidence bundle.

Prior downstream and benchmark evidence carries forward because every runtime
R implementation, DESCRIPTION and NAMESPACE are unchanged; the documentation
R source has identical parsed expressions. This is explicitly transferred
evidence, not a claim of fresh downstream checks or timing runs. The exact
rebuilt archive, checksums, CI source identities, and remaining NOTES are in the
release handoff supplied with this review.

Primary technical references:

- Erichson et al. (2020), [Sparse Principal Component Analysis via Variable Projection](https://arxiv.org/html/1804.00341v3), [DOI](https://doi.org/10.1137/18M1211350).
- R Core, [QR decomposition and its rank/pivot conventions](https://stat.ethz.ch/R-manual/R-devel/library/base/html/qr.html).
- R Core, [prcomp and preprocessing conventions](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html).
