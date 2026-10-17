# Independent review of objective-geometry tests

Reviewed `tests/testthat/test-objective-geometry.R` on 2026-10-17.
No implementation or repository test files were edited by this reviewer.
The bounded experiment executes only the first two blocks, saved in
`focused-objective-geometry.R`, against the installed checked rc2 package.

## Mathematical oracle assessment

For fixed B, maximizing the reconstruction cross-term over orthonormal A
has value the nuclear norm of X' X B. Therefore

`0.5 ||X||_F^2 + 0.5 ||X B||_F^2 - ||X' X B||_*`

is the correct reduced reconstruction objective. Central differences of
this scalar function supply a gradient independently of the production
thin gradient or its computed polar A. The explicit ridge derivative and
relative penalty coefficients then determine a complete one-step update
for either proximal penalty. This oracle directly covers the algorithm's
update, rather than merely checking that its final objective decreased.

The controlled starting point is safely inside the smooth domain: the
smallest singular value of X' X B is 21.14345. For the l0 comparison, the
nearest difference between inactive and active scalar proximal costs is
0.01112401, so the discrete active-set oracle is safely away from a tie.
The maximum coefficient discrepancy is 1.867817e-11 for both penalties.
The 1e-7 allowance has substantial backend-rounding headroom while
remaining far below the errors from the injected update faults.

The correlated-data KKT test uses the full residual and separate active
and zero-coefficient l1 conditions. Those conditions are sufficient for
the convex B subproblem at fixed A, but do not establish joint global
optimality. Checking the nuclear-norm maximum is an independent value
reference for the orthogonal Procrustes subproblem. The requirement that
both zero and nonzero coefficients occur prevents a vacuous subgradient
branch. The test correctly treats objective convergence tolerance and
coefficient/gradient tolerance as different quantities.

Two small improvements were recommended to the root author:

1. Use the literal requested alpha and beta multiplied by the independently
   calculated leading squared singular value in the KKT oracle. Using the
   returned effective coefficients is internally coherent but permits a
   consistent scaling error in the optimizer and reporting to survive that
   individual check. The first finite-difference block already catches it.
2. Add an explicit A' A = I assertion to the penalized correlated fit. The
   nuclear-norm gain equality characterizes the constrained optimum only
   among feasible A; the scalar equality alone does not prove feasibility.

The third block's change-of-units relation is correct. A negative global
multiplier leaves relative loading penalties unchanged, scales their
absolute coefficients and objective quadratically, scales scores and
robust errors linearly, and requires gamma to scale by its absolute value.
Reconstruction comparisons appropriately avoid demanding identical factor
signs. Shared RNG seeds are necessary for the randomized comparison.

## Isolated fault injections

`geometry-mutation-checks.R` temporarily replaces only `.spca_solve` inside
one R process, restoring it after every case. Installed files and
repository implementation files remain unchanged. Raw per-case testthat
objects and compact CSVs are retained.

| Fault | Finite-difference failures | KKT-block failures | Induced solver warnings |
|---|---:|---:|---:|
| Unmodified control | 0 of 2 | 0 of 5 | 0 |
| Omit ridge shrinkage | 2 of 2 | 3 of 5 | 1 |
| Reverse smooth-gradient sign | 2 of 2 | 3 of 5 | 1 |
| Treat alpha as absolute in proximal threshold | 2 of 2 | 3 of 5 | 0 |

There were no test errors or skips. Warnings in two injected KKT cases are
expected consequences of the intentionally incorrect optimizer's objective
increase; the unmodified test run has no warnings. Both independent blocks
catch each of the three mistakes, not merely the combined file.

These tests improve evidence for the implemented model substantially. They
cannot certify a global solution to the nonconvex joint problem, statistical
feature recovery, or all floating-point environments, and they do not make
such claims.

## Final test refinement and evidence transfer

The final KKT oracle uses independently calculated effective penalties and
one explicit orthogonality assertion. A parsed-expression comparison verified
that the preserved mutation snapshot (2 + 5 assertions) already included the
independent penalties and differs from the final block (2 + 6 assertions)
only by the added orthogonality assertion. The mutation experiment was not
rerun on the strengthened block. Existing mathematical conditions were
retained, the remaining oracle computation is unchanged,
and the complete final suite was executed separately. The mutation results
are therefore transferred evidence, not a fresh test of the final archive.
The versioned mutation scripts accept repository and output-directory
arguments; the original executed harnesses remain in the evidence bundle.
