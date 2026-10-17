# Independent second-pass solver review

Reviewed source: `erichson/spca`, merged commit
`62412f258465829d3f7110e1050d739374103781` (tree
`b0afe20f4a2e87f880600c1519fd9c5ea1b8f3e7`). Date: 2026-10-17.
The review reread the source, the existing tests, and the primary paper rather
than treating the earlier mathematical note or passing tests as proof.

Files in scope: `R/solver.R`, `R/spca.R`, `R/rspca.R`, `R/robspca.R`,
`R/utils.R`, with `R/preprocess.R` and `R/spca_path.R` read for context.
The only repository edit made by this reviewer is the new test file
`tests/testthat/test-robust-reference.R`.

## Conclusion

I found no additional confirmed mathematical defect in the fitting updates.
The Gram-equivalent reduction, polar-factor orientation, normalized penalty
scaling, elastic-net step, hard threshold, robust thin gradient, and sparse-error
update are mutually consistent with the stated objective. This conclusion is
supported by an explicit derivation and fresh independent scalar references;
it does not follow only from the existing tests.

Two qualifications matter: objective stagnation is not a stationarity guarantee,
and the projected value function is nonsmooth at rank loss. Neither invalidates
the actual fixed-block update implemented here. The strongest immediate
improvement is additional independent mathematical tests and better diagnostics,
not a new numerical backend or a broad feature expansion.

## 1. Deterministic and randomized update derivation

Let C be the right SVD factor D V' of the data or its sketch. Its Gram matrix
is unchanged: C'C equals X'X for the deterministic method and the sketch Gram
matrix for the randomized method. Therefore it preserves every term of the
penalized reconstruction objective, not just the PCA initialization.

For fixed B, expanding the squared residual and using A'A = I gives

    0.5 ||C - C B A'||F^2
      = 0.5 ||C||F^2 + 0.5 ||C B||F^2 - tr(A' C' C B).

The thin polar factor of C'C B maximizes the trace over A'A = I. Thus the
SVD at `solver.R:43` and `A = U V'` at line 44 have the correct orientation.
The same expansion justifies the fast objective at lines 4–19. Its explicit
blocked-residual fallback addresses cancellation near exact reconstruction;
a separate reviewer is probing that branch.

For fixed A, include the ridge term in the smooth loss. The gradient in B is

    g = C'C(B - A) + beta_eff B,

and its Lipschitz constant is

    L = ||C||2^2 + beta_eff = ||C||2^2 (1 + beta).

The actual update at lines 48–49 is exactly prox_{alpha_eff/L P}(B - g/L).
Here P(B) is the entrywise l1 norm or the number of nonzero coefficients.
The explicit normalization in the code avoids unnecessarily forming a large
scaled gradient; it does not change the penalized problem.

For l1, the proximal minimizer is soft thresholding at alpha_eff/L. For l0,
it is hard thresholding at sqrt(2 alpha_eff/L); selecting zero at equality
is a valid minimizer of the scalar proximal subproblem. The code uses these
thresholds with the ridge term kept in the smooth gradient. In particular,
the l0 branch is not an arbitrary cutoff applied after fitting.

The quadratic Taylor majorizer with this L dominates the fixed-A smooth loss.
The globally minimized separable proximal subproblem cannot increase that
majorizer, and the exact A update cannot increase the objective either. This
establishes blockwise nonincrease in exact arithmetic for both penalties,
without requiring the l0 penalty to be convex.

## 2. Robust updates and independent Huber oracle

Writing Y = X - S, the robust fixed-B,S A subproblem maximizes
tr(A'Y'XB). `solver.R:88` therefore uses the correct cross-product. The smooth
B gradient is

    X'(XB - (X-S)A) + beta_eff B.

It has the same Lipschitz bound as the ordinary fit because A'A = I. The thin
expression at lines 90–91 is algebraically correct and avoids constructing a
full n-by-p gradient residual. At fixed A,B, S = soft(X-XBA', gamma) is the
exact minimizer of the sparse-error subproblem. The objective must then use
X-XBA'-S; lines 93–98 do so.

A new independent test profiles S out analytically via the Huber loss. For
one variable, write b = |B| and choose A's sign optimally. No minimum requires
b outside [0,1], and the elastic-net scalar objective is

    sum_i H_gamma(x_i (1-b))
      + alpha ||x||2^2 b + 0.5 beta ||x||2^2 b^2.

Its derivative on [0,1] is

    -sum_i |x_i| min(|x_i| (1-b), gamma)
      + alpha ||x||2^2 + beta ||x||2^2 b.

A one-dimensional root solve gives an independent optimum with no reuse of the
package's alternating updates, initialization, or sparse-error implementation.
For x = (-5,-2,-1,0.5,1,3,20), alpha = 0.01, beta = 0.05, gamma = 0.5,
center = FALSE, the fitted coefficient is approximately 0.5295218, versus
0.5295209 from the root oracle. The fitted S is nonzero. At tol = 1e-13 the
coefficient error is about 9e-7; the unit-test coefficient bound is 5e-6 to
accommodate objective-based stopping and ordinary numerical-library variation.
The objective is checked much more closely against the scalar minimum.

For l0, the active scalar branch has the same derivative with alpha = 0.
Comparing its objective plus the l0 charge against b=0 supplies a complete
scalar global reference. The test checks a safely active case and a safely
inactive case whose default initialization reaches the corresponding basin.
It does not generalize that global-optimum result to multidimensional l0 fits.

For gamma=0, the error matrix can absorb all data for free. With positive
loading penalties the exact solution B=0,S=X has objective zero. The new test
checks this on nonzero data for both penalties. When gamma=alpha=beta=0, every
B paired with S=X-XBA' has zero objective; failure to select B=0 would not be
a mathematical error because the solution is nonunique.

Executed new test file: R 4.6.1, one BLAS/OpenMP thread, 2 test blocks,
23 assertions passed, no failures, warnings, or skips. Log:
`robust-reference-test-r461.log`. Additional scalar probes and raw output:
`scalar-robust-and-final-polar.R`, `scalar-robust-results.csv`.

## 3. Returned A and objective convergence

At the end of a cycle A was optimized for the preceding B (and preceding S
in the robust fit). The returned A,B,S are a valid feasible tuple, and the
recorded objective is the objective of that tuple. A is not necessarily the
exact conditional minimizer for the final B,S. This matters particularly if
max_iter is small or a loose tolerance stops early.

I measured the effect on ten independently seeded contaminated 60-by-8
problems, with k=2, alpha=0.02, beta=0.01 and gamma=0.5 for robust fits. At the
default tol=1e-5, recomputing only the final polar factor reduced the objective
by at most about 1.69e-6 relative for spca and 2.16e-6 for robspca. The maximum
normalized tangent-space A residuals were approximately 0.00116 and 0.00127.
At tol=1e-10, the relative loss gaps fell below 1.96e-11 and 3.77e-11, with
A residuals below 3.51e-6 and 5.35e-6. The raw per-seed observations are in
`final-polar-results.csv`.

This is consistent with objective changes becoming second-order near a
stationary point: a tolerance on function changes is not the same tolerance
on parameters or first-order residuals. The documentation already says that
converged means objective-based stopping. I do not recommend changing the
iteration ordering or silently adding an A/S polishing pass merely to make
A conditionally optimal on return. That would alter historical trajectories
and objective bookkeeping without fixing a demonstrated defect.

A useful extension would be an explicitly named first-order/proximal residual
at the returned tuple, perhaps behind a diagnostics option. It must include
all relevant blocks, distinguish the l0 fixed-point interpretation, use the
sketched objective for rspca, and avoid claiming global optimality.

## 4. Rank loss and the limits of the theoretical claim

Eliminating A yields the exact value function

    v(B) = 0.5 ||C||F^2 + 0.5 ||CB||F^2 - ||C'C B||*,

where ||.||* is the nuclear norm. Its polar factor is unique when C'C B has
full column rank, but not generally at rank loss. In particular v is
nondifferentiable at B=0 unless the data are zero. One cannot infer a globally
Lipschitz gradient of v from the fixed-A least-squares Lipschitz bound.

This does not undermine the source derivation above: the implementation uses
a fixed-A quadratic majorizer at each B step, for which the spectral bound
is valid even at rank loss. The exact SVD supplies an admissible polar factor
when it is nonunique. Nonconvexity and the choice among polar factors can
still affect the attained stationary point or support.

The primary paper's Algorithm 1 and Algorithm 2 provide the relevant
variable-projection and robust formulation. Its Appendix C.3 invokes
convexity of psi in its displayed proof; that argument should not be presented
as a blanket guarantee for every l0/rank-deficient case. The package should
promise correctly implemented updates and transparently defined stopping,
not global optimality or uniqueness. The existing documentation is broadly
consistent with this qualification.

Primary reference read for this review:
Erichson et al., Sparse Principal Component Analysis via Variable Projection,
SIAM Journal on Applied Mathematics 80(2), 977–1002 (2020),
[DOI 10.1137/18M1211350](https://doi.org/10.1137/18M1211350),
[author manuscript](https://arxiv.org/html/1804.00341v3), especially equations
(10)–(19), Algorithms 1–2 and Appendix C.3. The nonsmoothness observation above
also follows directly from the displayed nuclear-norm identity.

## 5. Prioritized recommendations

1. Keep the new active-Huber and gamma=0 tests. They resolve an actual gap:
   the previous active robust regression reconstructed the objective but did
   not independently locate the optimum; the previous analytical penalty
   references only covered an orthogonal non-robust regime.
2. Add an independent correlated-design loading KKT reference and carefully
   defined convergence diagnostics if feasible. Another reviewer is handling
   those paths. Avoid assertions that objective tolerance equals parameter
   tolerance or that arbitrary warm starts reproduce principal directions.
3. Retain the current SVD/Matrix infrastructure. No evidence from this review
   justifies replacing the backend, normalizing fitted B after optimization,
   or sparsifying the robust workspace immediately.
4. A cached Gram-matrix path for small p could reduce repeated products in
   deterministic fits, but costs another p-by-p allocation and squares the
   conditioning of an explicit covariance representation. It should be a
   separate measured, size-limited optimization, not an unbenchmarked default.
5. A fully sparse robust workspace would be a larger algorithm/memory project:
   S must be updated from a dense residual and the public result currently
   returns the full error matrix. It should wait for a clear workload and
   independent blocked/sparse equivalence and memory measurements.
6. Reconstruct/inverse-transform convenience or multistart selection can be
   useful follow-on APIs, but the core mathematics does not require them.
   Better convergence interpretation and independent references have higher
   priority for this release than broad feature expansion.
