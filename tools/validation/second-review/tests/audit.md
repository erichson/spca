# Second review: unit-test adequacy

Date: 2026-10-17. Review baseline: merged commit
`62412f258465829d3f7110e1050d739374103781` (sparsepca 1.0.0).
All five original test files and the harness were read alongside every R
implementation file. The review distinguishes mathematical evidence from
implementation agreement and targets important gaps instead of maximizing
assertion count.

## Assessment

The existing 29 blocks are a strong foundation. The robust regression
reconstructs the intended objective from returned factors; the PCA checks
match preprocessing; the randomized checks use spectral-tail bounds and
subspaces; the penalty checks use analytic scalar optima. None of those
should be discarded. The former uninvoked/signed-residual tests have already
been replaced.

Remaining weaknesses identified:

1. Sparse/dense fit comparisons use p=12, k=3, o=20: the sketch is full width.
   They establish storage/preprocessing equivalence, but do not establish
   compressed sparse-sketch quality. The q regression exercises only the
   dense rsvd backend. Assigned to the independent sparse-operator review.
2. Scalar penalty references are orthogonal. They catch threshold scaling,
   but do not strongly exercise a correlated iterative gradient. A reduced
   objective finite-difference check and fixed-factor KKT reference are
   assigned to the root/solver review.
3. Sparse preprocessing references call the same package helper on dense
   and sparse input. This is useful cross-branch agreement, but shared
   preprocessing errors can survive. Independent base scale/sum-of-squares
   references complement them; the existing deterministic prcomp test
   already supplies such an oracle for ordinary finite dense inputs.
4. The adjusted-variance oracle uses qr.resid, a different calculation but
   the same QR family as production. A closed-form geometry example is
   valuable particularly when zero/dependent columns trigger reordering.
5. No targeted test forced cancellation in the optimized residual identity
   or explicitly distinguished a significant objective increase from
   numerical stagnation. These are numerical safeguards, so success on
   benign end-to-end examples is insufficient evidence.

## Existing coverage map

| Test block | Risk or contract checked | Scope of its evidence |
|---|---|---|
| Implicit centered/scaled products | Forward and adjoint values; structural Matrix coercion; adjoint identity | Dense branch shares preprocessing; many storage classes are covered |
| Fully stored columns with offsets | Cancellation in centered sparse products | Excellent adversarial offset case; helper level |
| Sparse/dense fits | Scores, reconstruction, total/adjusted energy | Full-width randomized sketch, not compressed quality |
| Sparse allocation profile | Full n by p dense intermediate | Checks large individual allocations, not total or peak process memory |
| Sparse missing/unit diagonal | Row omission and implicit structural entries | Correctly uses residual energy for tied bases |
| Scalar elastic net | Penalized magnitude and effective coefficients | Independent analytic optimum; orthogonal design |
| Scalar l0 | Active/inactive global scalar optima and true objective | One-component conditional reference |
| Default elastic net | Default argument compatibility | Intentionally a same-implementation comparison |
| Zero component variance | Component order after leading zero | Good targeted regression; extended below for additional dependence |
| Loading warm start | Reused objective and score convention; sign ambiguity | Tests numerical meaning, not arbitrary matching signs |
| Invalid penalty/start | Unsupported controls and shapes | Broad errors protect public failure contract |
| Unpenalized PCA | Reconstruction, orthogonality, sdev, scores and total energy | Independent prcomp plus explicit preprocessing, tall/wide |
| Randomized exact/tail | Exact rank and spectral-tail quality | Two seeds per controlled case; appropriate relative bound |
| Repeated spectrum | Invariant subspace and optimal residual | Avoids nonunique-vector expectations |
| Rank deficiency | Output shapes, finite diagnostics, residual | Covers all three solvers |
| Correlated variance | Sequential adjusted versus marginal energy | qr.resid is partially independent, not a separate linear algebra backend |
| Ridge magnitude | Retained loading shrinkage and variance | Independent closed form |
| Robust objective | Correct R-S residual, soft threshold, descent, iteration count | Strong regression exposes old premature stop |
| Large gamma | Robust/ordinary limiting equivalence | Tests no-outlier regime |
| q effect | Honored power-iteration argument and quality gain | Dense backend, one deliberate regression seed |
| Zero/constant input | Finite early return, shape, summary convention | Does not exercise general zero objective stopping |
| Single variable | Dimension retention, names, variance | Covers all solvers |
| Missing rows | Warning, names, fitted means | Covers ordinary numeric inputs |
| Invalid numerical inputs | Numeric/scalar/domain validation | Does not prove rejection precedes every expensive preprocessing operation |
| Iteration/rank bounds | max_iter=1 nonconvergence and rank clipping | Preserves deliberate public behavior |
| Prediction preprocessing | Stored centering/scaling, reorder, Matrix, one/zero rows | Broad useful output reference |
| Prediction errors | Named-feature mismatch and missing new values | Covers common misuse |
| Deterministic path | Alpha order, warm starts, counts, print | Numerical equivalence alone cannot prove one factorization; source review does |
| Randomized path/RNG | One sketch stream; deterministic methods leave seed intact | Direct RNG-state evidence complements source review |

## Targeted additions made in this review

Only `tests/testthat/test-numerical-safeguards.R` was edited by this reviewer.
It adds three blocks, 19 assertions:

* A general rank-two reconstruction and a near-exact full reconstruction
  compare `.spca_loss` to the actual represented residual. In the latter
  case the true loss is approximately 4.25e-17, while subtracting the three
  large energy terms gives zero. A ratio check prevents an absolute
  tolerance from making the test vacuous. The relative tolerance is 1e-6
  to allow harmless rounding changes in stable residual implementations.
* Stopping tests distinguish meaningful descent, a roundoff-scale increase,
  a significant increase, and a zero-to-zero objective on nonzero-energy
  data. The existing zero-data shortcut does not cover that denominator.
* A known orthogonal input and loading geometry gives adjusted score
  energies (0,2,0,8,18) and marginal energies (0,2,18,10,22), before division
  by n-1. The production QR pivots in order (2,4,5,1,3); the expected values
  come directly from orthogonal axes, without calling any QR oracle.

These are helper-level unit tests because the risk is specifically in the
shared numerical safeguards and result formatter. End-to-end solver tests
remain responsible for proving those helpers are connected correctly.

## Executed validation and fault injection

R 4.6.1, single-threaded OpenBLAS, installed checked rc2 library. The new
file passes all 19 assertions with no failures, errors, warnings, or skips.
`unmodified-results.rds/csv` and `sessionInfo.txt` preserve the results.
The initial CSV export attempted to write testthat's list-column result;
it failed after the tests had passed. The export was corrected by removing
that result column; the full object remains in the RDS. This was an evidence
serialization error, not a test failure.

Four isolated process-local namespace mutations establish that the new
checks actually detect the intended faults. Installed package files and
repository source were never modified by these injections.

| Injected fault | Passed assertions | Failed assertions | Errors |
|---|---:|---:|---:|
| Unmodified control | 19 | 0 | 0 |
| Remove direct-residual fallback | 18 | 1 | 0 |
| Mark significant objective increase converged | 18 | 1 | 0 |
| Ignore data energy in roundoff allowance | 17 | 2 | 0 |
| Forget mapping QR diagonal energies through the pivot | 18 | 1 | 0 |

All runs have zero unexpected warnings and skips. These targeted mutation
checks are development evidence, not additions to the CRAN test suite.
`mutation-checks.R`, `mutation-summary.csv`, and per-case RDS/CSV outputs
make the experiment reproducible.

## Recommendations deliberately deferred

Do not add timing assertions, expectations of identical randomized factors,
full duplicate reference solvers, or blanket parameter grids. They would
increase maintenance cost without resolving a named uncertainty. Broader
hard-threshold global-optimality claims are inappropriate for this nonconvex
algorithm. Exact cardinality, automatic statistical tuning and cross-validation
need separate API/design work rather than opportunistic unit-test additions.
