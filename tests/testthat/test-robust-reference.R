test_that("active robust fits agree with an independently profiled Huber optimum", {
  # For one variable A is +1 or -1. Minimizing over A makes b = |B|
  # nonnegative, and no optimum needs b > 1. Eliminating S gives a convex
  # scalar Huber problem for the elastic-net case. This reference does not
  # reuse the package's alternating updates or sparse-error calculation.
  x <- c(-5, -2, -1, 0.5, 1, 3, 20)
  X <- matrix(x, ncol = 1)
  gamma <- 0.5
  beta <- 0.05
  spectral2 <- sum(x^2)
  huber <- function(residual) {
    a <- abs(residual)
    sum(ifelse(a <= gamma, 0.5 * a^2, gamma * a - 0.5 * gamma^2))
  }
  scalar_minimum <- function(alpha) {
    derivative <- function(b) {
      -sum(abs(x) * pmin(abs(x) * (1 - b), gamma)) +
        alpha * spectral2 + beta * spectral2 * b
    }
    if (derivative(0) >= 0) return(0)
    stats::uniroot(derivative, c(0, 1), tol = 1e-13)$root
  }

  alpha <- 0.01
  expected_b <- scalar_minimum(alpha)
  expected_objective <- huber(x * (1 - expected_b)) +
    alpha * spectral2 * expected_b + 0.5 * beta * spectral2 * expected_b^2
  fit <- robspca(X, 1, alpha = alpha, beta = beta, gamma = gamma,
                center = FALSE, tol = 1e-13, max_iter = 2000, verbose = FALSE)
  expect_true(fit$converged)
  expect_true(any(fit$sparse != 0))
  # Objective-based stopping controls function values, not coefficient error.
  # The measured scalar error is about 1e-6; this bound allows ordinary BLAS
  # rounding variation while excluding an incorrectly scaled/updated solution.
  expect_lt(abs(abs(fit$loadings[1, 1]) - expected_b), 5e-6)
  expect_equal(tail(fit$objective, 1), expected_objective, tolerance = 1e-9)
  expect_equal(huber(as.vector(X - tcrossprod(fit$scores, fit$transform))) +
                 fit$alpha * sum(abs(fit$loadings)) +
                 0.5 * fit$beta * sum(fit$loadings^2),
               tail(fit$objective, 1), tolerance = 1e-12)

  # With l0 regularization the active optimum has no alpha derivative.
  # Compare its value with b = 0 to select the global scalar minimum. These
  # two cases deliberately start in the corresponding attraction region;
  # generic multidimensional l0 fits are not promised to find a global optimum.
  active_b <- scalar_minimum(0)
  active_base <- huber(x * (1 - active_b)) +
    0.5 * beta * spectral2 * active_b^2
  zero_value <- huber(x)
  for (alpha in c(0.001, 0.5)) {
    active_value <- active_base + alpha * spectral2
    expected_b <- if (active_value < zero_value) active_b else 0
    expected_objective <- min(active_value, zero_value)
    fit <- robspca(X, 1, alpha = alpha, beta = beta, gamma = gamma,
                  center = FALSE, penalty = "l0", tol = 1e-13,
                  max_iter = 2000, verbose = FALSE)
    expect_true(fit$converged)
    expect_true(any(fit$sparse != 0))
    expect_lt(abs(abs(fit$loadings[1, 1]) - expected_b), 5e-6)
    expect_equal(tail(fit$objective, 1), expected_objective, tolerance = 1e-9)
  }
})

test_that("zero robust error penalty allows the sparse term to absorb all data", {
  X <- matrix(c(-5, -2, -1, 0.5, 1, 3, 20), ncol = 1)
  # gamma = 0 removes the cost of S. Positive loading penalties then have
  # the exact minimum B = 0, S = X, with objective zero for either penalty.
  for (penalty in c("elasticnet", "l0")) {
    fit <- robspca(X, 1, alpha = 0.5, beta = 0.1, gamma = 0,
                  center = FALSE, penalty = penalty, max_iter = 100,
                  tol = 1e-12, verbose = FALSE)
    expect_true(fit$converged)
    expect_equal(fit$loadings, matrix(0, 1, 1), ignore_attr = TRUE)
    expect_equal(fit$scores, matrix(0, nrow(X), 1), ignore_attr = TRUE)
    expect_equal(fit$sparse, X, ignore_attr = TRUE)
    expect_equal(tail(fit$objective, 1), 0)
  }
})
