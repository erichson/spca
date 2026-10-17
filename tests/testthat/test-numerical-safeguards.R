test_that("residual loss remains accurate near an exact reconstruction", {
  C <- diag(c(8, 4, 2, 1))
  cases <- list(
    list(A = cbind(rep(0.5, 4), c(0.5, -0.5, 0.5, -0.5)),
         B = matrix(c(0.5, 0, -0.2, 0.1, 0, 0.3, 0.1, -0.4), 4, 2)),
    list(A = diag(4), B = (1 - 1e-9) * diag(4))
  )
  for (case in cases) {
    A <- case$A
    B <- case$B
    CA <- C %*% A
    CB <- C %*% B
    # Evaluate the represented residual directly. At the near-exact fit,
    # subtracting its three large energies rounds the loss to zero.
    expected <- 0.5 * sum((C - CB %*% t(A))^2)
    observed <- sparsepca:::.spca_loss(C, CA, CB, A, sum(C^2))
    expect_gt(expected, 0)
    # Compare relatively: an absolute tolerance would let a spurious zero
    # pass when the true loss is approximately 4.25e-17. One part per million
    # allows harmless rounding differences in stable residual evaluations.
    expect_equal(observed / expected, 1, tolerance = 1e-6)
  }
})

test_that("stopping distinguishes descent, roundoff, and objective increases", {
  progress <- sparsepca:::.spca_progress
  decreasing <- progress(2, 1, tol = 0.1, evaluation_scale = 1000)
  expect_equal(decreasing$improvement, 1)
  expect_false(decreasing$stop)
  expect_false(decreasing$converged)

  # This increase is below the rounding uncertainty of an objective
  # evaluated from a data matrix with squared norm 1000.
  expect_warning(tiny <- progress(1, 1 + 1e-12, tol = 0,
                                  evaluation_scale = 1000), NA)
  expect_true(tiny$stop)
  expect_true(tiny$converged)

  expect_warning(increased <- progress(1, 1.001, tol = 1e-5,
                                       evaluation_scale = 1000),
                 "objective increased")
  expect_true(increased$stop)
  expect_false(increased$converged)
  expect_lt(increased$improvement, 0)

  # Nonzero-data fits can also attain zero objective; the zero-data fast
  # return alone does not exercise this denominator.
  zero <- progress(0, 0, tol = 0, evaluation_scale = 1000)
  expect_identical(zero$improvement, 0)
  expect_true(zero$stop)
  expect_true(zero$converged)
})

test_that("adjusted variances retain geometry after zero and dependent scores", {
  # The five input axes are centered, mutually orthogonal, and have
  # squared norm 2. This gives a closed-form oracle without another QR.
  X <- rbind(diag(5), -diag(5))
  B <- cbind(c(0, 0, 0, 0, 0), c(1, 0, 0, 0, 0),
             c(3, 0, 0, 0, 0), c(1, 2, 0, 0, 0), c(1, 1, 3, 0, 0))
  prepared <- sparsepca:::.spca_preprocess(X, center = FALSE, scale = FALSE)
  fitted <- list(B = B, A = diag(5), objective = 0,
                 converged = TRUE, improvement = 0)
  result <- sparsepca:::.spca_result(fitted, prepared, "spca", 0, 0,
                                    "elasticnet")
  # The third score duplicates the second. The fourth adds 2*axis2 and
  # the fifth adds 3*axis3; these contributions stay in component order.
  expect_equal(result$eigenvalues * (nrow(X) - 1), c(0, 2, 0, 8, 18),
               tolerance = 1e-12)
  expect_equal(result$score_variance * (nrow(X) - 1), c(0, 2, 18, 10, 22),
               tolerance = 1e-12)
})
