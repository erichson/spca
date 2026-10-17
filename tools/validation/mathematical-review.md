# Independent numerical review for sparsepca 1.0.0

## Scope and baseline

The review compared all three original R implementations with the variable
projection formulation in Erichson et al., *SIAM Journal on Applied Mathematics*
80(2), 977–1002 (2020), DOI
[10.1137/18M1211350](https://doi.org/10.1137/18M1211350), available as
[arXiv:1804.00341](https://arxiv.org/abs/1804.00341). The published 0.1.2 R sources,
GitHub revision `fa896289a68801a4798f53948ff6c2b8f7c54e1c`, and the R sources in
the maintainer's attached archive were identical. Baseline numerical probes ran
under R 4.6.1 on Ubuntu 24.04.3 with one OpenBLAS thread.

## Confirmed corrections

### Robust residual and premature stopping

The robust criterion is

`0.5 * ||X - X B A' - S||_F^2 + alpha_eff * ||B||_1 +
0.5 * beta_eff * ||B||_F^2 + gamma * ||S||_1`.

If `R = X - X B A'`, the squared residual is `R - S`. The original expression
charged the sparse-outlier penalty while leaving the outliers in the quadratic
loss. The corrected expression agrees with equations 17–19 of the paper.

The 80 × 8 contaminated rank-two example now stored in the regression test
stopped after two iterations in the published implementation. Its recorded
objective was approximately 190.5643, whereas evaluation of the intended
criterion at that same returned fit gave 43.37158. Applying only the reported
one-line correction produced 1,933 iterations at `tol = 1e-10`, with a final
objective of approximately 31.51068 and a monotonically decreasing trace.
The regression suite checks the true objective, the soft-thresholded outlier
matrix, continued progress beyond two iterations, and monotonicity.

### Randomized power iterations

The original `rspca` accepted `q` but always called `rqb(..., q = 2)`. Fixed-seed
fits with `q = 0` and `q = 4` therefore produced exactly the same reconstruction.
The regression uses a controlled decaying spectrum and checks that additional
power iterations change the fit and reduce approximation error.

### Variance diagnostics

The previous `eigenvalues` were the singular values of a Procrustes
cross-product from before the final loading update, divided by `n - 1`.
Those quantities are not variances of the returned scores. For a centered
orthogonal example with leading squared singular value 32 and `beta = 0.5`,
the correct loading norm is 2/3. The old reported value was 4.266667 while the
actual score variance was 2.844444.

The release uses sequential QR-adjusted score energies divided by `n - 1`,
with numerically dependent columns mapped back to their original component
positions. This adjustment follows the motivation in Section 3.4 of
[Zou, Hastie and Tibshirani (2006)](https://hastie.su.domains/Papers/spc_jcgs.pdf).
The historic field name `eigenvalues` is retained, but the documentation makes
clear that these are not generally covariance eigenvalues. Unadjusted
per-score energies are separately available as `score_variance`.

When centering is disabled, both the numerator and total `var` use uncentered
second moments with the same `n - 1` divisor. This preserves agreement with
`prcomp(center = FALSE)` at zero penalties and fixes the old inconsistent
variance ratios. The tests compare both centered and uncentered cases with
matching preprocessing, and independently compute the sequential adjustment
by regressing each score on the preceding scores.

### Degeneracies and termination

The old zero-matrix calculation eventually evaluated `0 / 0` in its stopping
rule; `max_iter = 0` accessed an undefined SVD update. The release handles zero
energy explicitly and validates iteration counts, ranks, finite real inputs,
and tuning parameters. The objective-increase tolerance uses the data-energy
scale of its floating-point evaluation, so machine-level residual changes near
exact fits are not falsely reported as failed convergence. Meaningful increases
are still distinguished from convergence.

Rank-one dimension dropping was initially identified as a risk, but a real R
probe showed that the original one-variable case worked because subsequent
matrix products restored its dimensions. It is not listed as a confirmed
historical bug. Explicit matrix shapes are now tested, including that case.

## Solver equivalence and memory reduction

Let `C = D V'` be the right factor of an SVD of the preprocessed data or its
randomized sketch. Then `C' C = X' X` for the deterministic solver. The
Procrustes update is the polar factor of `C' C B`; the negative smooth gradient
is `C' C (A - B) - beta_eff B`. With
`alpha_eff = alpha * ||C||_2^2` and
`beta_eff = beta * ||C||_2^2`, the implemented elastic-net step

`soft((B + C' C (A - B) / ||C||_2^2) / (1 + beta), alpha / (1 + beta))`

is algebraically the original proximal-gradient step. The code retains the
estimated coefficient magnitudes; it does not normalize the loading columns
after optimization.

Because `A' A = I`, the quadratic residual can be evaluated as

`||C||_F^2 - 2 * <C A, C B> + ||C B||_F^2`.

This removes the large residual from ordinary iterations. Near an exact
reconstruction, subtraction of large energies is inaccurate; the implementation
then evaluates the actual residual in bounded row blocks. The deterministic
initialization also avoids calculating unused left singular vectors.

For robust fitting, the equivalent thin gradient uses `X' ((X - S) A - X B)`.
The dense outlier matrix is still required by the public returned object and
the robust objective. This release does not claim that robust SPCA is a
matrix-free sparse-input algorithm.

## l0 extension

The optional `l0` penalty uses the number of nonzero coefficients plus the
existing ridge term. The ridge term remains in the smooth gradient, and the
proximal threshold is `sqrt(2 * step * alpha_eff)`. Thus the new branch solves
the stated zero-norm-penalized subproblem; it does not apply an arbitrary hard
cutoff to a completed PCA fit. Tests use orthogonal scalar problems with known
active and inactive optima, and check objective descent. The objective is
nonconvex, so neither branch promises a global optimum.

## Validation boundaries

The mathematical suite checks true residual norms, matched PCA references,
known penalized optima, repeated-spectrum subspaces, approximately low-rank
spectral tails, zero and rank-deficient data, warm starts, and returned
diagnostics. It avoids signed sums of reconstruction errors and fixed signs
for singular vectors. Sparse operator, structural Matrix-class, prediction,
minimum-version, downstream, archive, and platform checks are recorded in
the broader maintainer and release reports. This note does not substitute
source review for those executed checks.
