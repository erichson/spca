test_that("the penalized update agrees with a finite-difference value-function gradient", {
  set.seed(802)
  X <- matrix(rnorm(18 * 6), 18, 6)
  X[, 2:6] <- X[, 2:6] + 0.4 * X[, 1:5]
  B <- matrix(rnorm(12), 6, 2) / 3
  # Eliminate A through the optimal Procrustes value (the nuclear norm).
  # This oracle never uses the implementation's thin gradient or polar A.
  # The chosen B gives a full-column-rank cross-product, away from the
  # nondifferentiability of the nuclear norm at rank loss.
  value <- function(B) {
    XB <- X %*% B
    0.5 * sum(X^2) + 0.5 * sum(XB^2) -
      sum(svd(crossprod(X, XB), nu = 0, nv = 0)$d)
  }
  h <- 1e-5
  gradient <- B
  for (j in seq_along(B)) {
    plus <- minus <- B
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    gradient[j] <- (value(plus) - value(minus)) / (2 * h)
  }
  spectral2 <- svd(X, nu = 0, nv = 0)$d[1]^2
  alpha <- 0.035
  beta <- 0.2
  step <- 1 / (spectral2 * (1 + beta))
  z <- B - step * (gradient + beta * spectral2 * B)
  for (penalty in c("elasticnet", "l0")) {
    expected <- if (penalty == "elasticnet") {
      sign(z) * pmax(abs(z) - step * alpha * spectral2, 0)
    } else {
      # Compare the scalar proximal costs at zero and at the active minimum.
      z * (0.5 * z^2 > step * alpha * spectral2)
    }
    fit <- spca(X, 2, alpha = alpha, beta = beta, center = FALSE,
                penalty = penalty, B_init = B, max_iter = 1, verbose = FALSE)
    # Central differences are approximate. The observed coefficient error
    # is about 2e-11; 1e-7 permits backend rounding without hiding a wrong step.
    expect_lt(max(abs(fit$loadings - expected)), 1e-7)
  }
})

test_that("a correlated penalized fit satisfies independent optimality conditions", {
  set.seed(1701)
  X <- matrix(rnorm(55 * 7), 55, 7)
  X[, 2:7] <- X[, 2:7] + 0.65 * X[, 1:6]
  fit <- spca(X, 3, alpha = 0.045, beta = 0.12,
              max_iter = 2000, tol = 1e-13, verbose = FALSE)
  D <- sweep(X, 2, fit$center, "-")
  A <- fit$transform
  B <- fit$loadings
  # Differentiate the full represented residual, not a compressed Gram
  # calculation. Active coefficients have a signed l1 derivative; at zero
  # the smooth derivative must lie in the l1 subgradient interval.
  spectral2 <- svd(D, nu = 0, nv = 0)$d[1]^2
  alpha_effective <- 0.045 * spectral2
  beta_effective <- 0.12 * spectral2
  gradient <- crossprod(D, (D %*% B) %*% t(A) - D) %*% A + beta_effective * B
  active <- B != 0
  residual <- gradient
  residual[active] <- gradient[active] + alpha_effective * sign(B[active])
  residual[!active] <- pmax(abs(gradient[!active]) - alpha_effective, 0)
  expect_true(fit$converged)
  expect_equal(crossprod(A), diag(ncol(A)), tolerance = 1e-12)
  expect_gt(sum(active), 0)
  expect_gt(sum(!active), 0)
  # Objective tolerance is not a coefficient/gradient tolerance. The
  # measured normalized KKT residual is about 1e-7 at tol = 1e-13.
  expect_lt(max(abs(residual)) / spectral2, 2e-6)

  cross <- crossprod(D, D %*% B)
  best_gain <- sum(svd(cross, nu = 0, nv = 0)$d)
  # The nuclear norm is the exact maximum over all orthonormal A. This
  # checks Procrustes optimality without comparing a particular SVD basis.
  expect_lt(abs(best_gain - sum(A * cross)) / sum(D^2), 1e-10)
})

test_that("relative loading penalties respect a change in data units", {
  set.seed(306)
  X <- matrix(rnorm(30 * 7), 30, 7)
  X[, 2:7] <- X[, 2:7] + 0.5 * X[, 1:6]
  X[4, 3] <- X[4, 3] + 9
  multiplier <- -3
  for (method in c("spca", "rspca", "robspca")) {
    for (penalty in c("elasticnet", "l0")) {
      args <- list(X = X, k = 2, alpha = 0.018, beta = 0.08,
                    center = TRUE, scale = FALSE, max_iter = 600,
                    tol = 1e-11, penalty = penalty, verbose = FALSE)
      if (method == "rspca") {
        args$o <- 1
        args$q <- 2
      }
      if (method == "robspca") args$gamma <- 0.6
      set.seed(121)
      original <- do.call(get(method), args)
      args$X <- multiplier * X
      # S changes units with X, so gamma must change by |multiplier|.
      # alpha and beta are relative to d1^2 and remain unchanged.
      if (method == "robspca") args$gamma <- abs(multiplier) * args$gamma
      set.seed(121)
      changed <- do.call(get(method), args)
      expect_true(original$converged && changed$converged)
      expected <- multiplier * tcrossprod(original$scores, original$transform)
      actual <- tcrossprod(changed$scores, changed$transform)
      expect_lt(norm(actual - expected, "F") / norm(expected, "F"), 1e-7)
      expect_equal(c(changed$alpha, changed$beta),
                   multiplier^2 * c(original$alpha, original$beta),
                   tolerance = 1e-12)
      expect_equal(tail(changed$objective, 1),
                   multiplier^2 * tail(original$objective, 1),
                   tolerance = 1e-9)
      if (method == "robspca") {
        expect_equal(changed$sparse, multiplier * original$sparse,
                     tolerance = 1e-7)
      }
    }
  }
})
