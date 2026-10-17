test_that("the robust objective includes the fitted sparse residual", {
  # This example stopped after two iterations in sparsepca 0.1.2. Merely
  # checking that an objective was returned would not detect the error.
  set.seed(1)
  u <- matrix(rnorm(80 * 2), 80, 2)
  v <- matrix(c(rep(1, 4), rep(0, 8), rep(1, 4)), 8, 2)
  X <- u %*% t(v) + matrix(rnorm(80 * 8, sd = 0.1), 80, 8)
  where <- cbind(c(3, 17, 29, 44, 62), c(1, 3, 4, 6, 8))
  X[where] <- X[where] + c(8, -10, 9, -8, 12)
  fit <- robspca(X, k = 2, alpha = 0.001, beta = 0.001, gamma = 0.5,
                max_iter = 200, tol = 1e-10, verbose = FALSE)
  xc <- sweep(X, 2, fit$center, "-")
  spectral2 <- svd(xc, nu = 0, nv = 0)$d[1]^2
  residual <- xc - fit$scores %*% t(fit$transform)
  expected <- 0.5 * sum((residual - fit$sparse)^2) +
    0.001 * spectral2 * sum(abs(fit$loadings)) +
    0.0005 * spectral2 * sum(fit$loadings^2) + 0.5 * sum(abs(fit$sparse))

  expect_true(any(fit$sparse != 0))
  expect_gt(length(fit$objective), 10)
  expect_equal(tail(fit$objective, 1), expected, tolerance = 1e-10)
  expect_true(all(diff(fit$objective) <= 1e-10 * max(fit$objective)))
  expect_lt(tail(fit$objective, 1), fit$objective[1] - 1)
  expect_equal(fit$sparse, sign(residual) * pmax(abs(residual) - 0.5, 0),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(crossprod(fit$transform), diag(2), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_identical(fit$iterations, as.integer(length(fit$objective)))
})

test_that("large gamma reduces robust SPCA to ordinary SPCA", {
  set.seed(118)
  X <- matrix(rnorm(40 * 8), 40, 8)
  ordinary <- spca(X, 3, alpha = 0.01, beta = 0.02,
                   max_iter = 1000, tol = 1e-10, verbose = FALSE)
  robust <- robspca(X, 3, alpha = 0.01, beta = 0.02, gamma = 1e8,
                   max_iter = 1000, tol = 1e-10, verbose = FALSE)
  expect_equal(robust$sparse, matrix(0, 40, 8), ignore_attr = TRUE)
  expect_equal(robust$scores %*% t(robust$transform),
               ordinary$scores %*% t(ordinary$transform),
               tolerance = 1e-7, ignore_attr = TRUE)
  expect_equal(tail(robust$objective, 1), tail(ordinary$objective, 1),
               tolerance = 1e-8)
})

test_that("q changes the randomized sketch and improves a difficult spectrum", {
  set.seed(701)
  left <- qr.Q(qr(matrix(rnorm(60 * 40), 60, 40)))
  right <- qr.Q(qr(matrix(rnorm(40 * 40), 40, 40)))
  spectrum <- 20 * seq_len(40)^(-0.7)
  X <- sweep(left, 2, spectrum, "*") %*% t(right)
  set.seed(22)
  quick <- rspca(X, 4, alpha = 0, beta = 0, center = FALSE,
                 o = 0, q = 0, verbose = FALSE)
  set.seed(22)
  refined <- rspca(X, 4, alpha = 0, beta = 0, center = FALSE,
                   o = 0, q = 4, verbose = FALSE)
  quick_reconstruction <- quick$scores %*% t(quick$transform)
  refined_reconstruction <- refined$scores %*% t(refined$transform)
  expect_gt(norm(quick_reconstruction - refined_reconstruction, "F"), 1e-4)
  expect_lt(sum((X - refined_reconstruction)^2),
            sum((X - quick_reconstruction)^2))
})

test_that("zero and constant data have finite, dimensionally valid fits", {
  for (fun in list(spca, rspca, robspca)) {
    for (X in list(matrix(0, 7, 4), matrix(rep(1:4, each = 7), 7, 4))) {
      fit <- fun(X, 2, scale = TRUE, verbose = FALSE)
      expect_identical(dim(fit$loadings), c(4L, 2L))
      expect_identical(dim(fit$transform), c(4L, 2L))
      expect_identical(dim(fit$scores), c(7L, 2L))
      expect_true(all(is.finite(fit$loadings)))
      expect_true(all(is.finite(fit$transform)))
      expect_equal(fit$scores, matrix(0, 7, 2), ignore_attr = TRUE)
      expect_equal(fit$eigenvalues, c(0, 0), ignore_attr = TRUE)
      expect_true(all(is.finite(fit$objective)))
      expect_true(all(is.finite(summary(fit))))
      expect_true(fit$converged)
    }
  }
})

test_that("one component and one variable keep their matrix dimensions", {
  X <- matrix(seq_len(12), ncol = 1, dimnames = list(letters[1:12], "value"))
  for (fun in list(spca, rspca, robspca)) {
    fit <- fun(X, 1, alpha = 0, beta = 0, verbose = FALSE)
    expect_identical(dim(fit$loadings), c(1L, 1L))
    expect_identical(dim(fit$transform), c(1L, 1L))
    expect_identical(dim(fit$scores), c(12L, 1L))
    expect_identical(rownames(fit$scores), rownames(X))
    expect_equal(fit$eigenvalues, var(as.vector(X)), ignore_attr = TRUE)
    expect_identical(dim(summary(fit)), c(4L, 1L))
  }
})

test_that("missing rows are omitted explicitly and preserve observation names", {
  X <- matrix(seq_len(24), 6, 4, dimnames = list(letters[1:6], LETTERS[1:4]))
  X[2, 3] <- NA_real_
  for (fun in list(spca, rspca, robspca)) {
    expect_warning(fit <- fun(X, 2, verbose = FALSE), "[Mm]issing|omit")
    expect_identical(nrow(fit$scores), 5L)
    expect_identical(rownames(fit$scores), rownames(X)[-2])
    expect_equal(fit$center, colMeans(X[-2, , drop = FALSE]))
  }
})

test_that("invalid numerical arguments fail before the optimizer", {
  X <- matrix(seq_len(24), 6, 4)
  invalid <- list(
    list(k = 0), list(k = -1), list(k = 1.5), list(k = NA_real_),
    list(k = Inf), list(k = c(1, 2)), list(alpha = -1), list(alpha = NA_real_),
    list(alpha = c(0.1, 0.2)), list(beta = -1), list(beta = Inf),
    list(max_iter = 0), list(max_iter = 1.5), list(tol = -1), list(tol = NA_real_),
    list(center = NA), list(scale = c(TRUE, FALSE)), list(verbose = 2)
  )
  for (fun in list(spca, rspca, robspca)) {
    for (args in invalid) {
      call_args <- utils::modifyList(list(X = X, k = 2, verbose = FALSE), args)
      expect_error(do.call(fun, call_args))
    }
    for (bad in list(matrix(1, 1, 4), matrix(numeric(), 0, 4),
                     matrix(numeric(), 4, 0), X + 1i,
                     matrix(Inf, 6, 4), matrix("a", 6, 4))) {
      expect_error(fun(bad, k = 1, verbose = FALSE))
    }
  }
  for (q in list(-1, 0.5, Inf, NA_real_)) {
    expect_error(rspca(X, 2, q = q, verbose = FALSE))
  }
  for (o in list(-1, 0.5, Inf, NA_real_)) {
    expect_error(rspca(X, 2, o = o, verbose = FALSE))
  }
  for (gamma in list(-1, Inf, NA_real_, c(1, 2))) {
    expect_error(robspca(X, 2, gamma = gamma, verbose = FALSE))
  }
})

test_that("iteration caps and legacy rank clipping have explicit diagnostics", {
  set.seed(191)
  X <- matrix(rnorm(30 * 6), 30, 6)
  for (fun in list(spca, rspca, robspca)) {
    fit <- fun(X, 2, alpha = 0.02, beta = 0.05, max_iter = 1, verbose = FALSE)
    expect_identical(fit$iterations, 1L)
    expect_length(fit$objective, 1)
    expect_false(fit$converged)
    clipped <- fun(X, 20, max_iter = 1, verbose = FALSE)
    expect_identical(dim(clipped$loadings), c(6L, 6L))
  }
})
