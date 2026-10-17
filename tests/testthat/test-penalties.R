test_that("elastic-net solutions agree with an orthogonal scalar reference", {
  X <- rbind(diag(c(4, 3, 2)), -diag(c(4, 3, 2)))
  alpha <- 0.05
  beta <- 0.1
  spectral2 <- 32
  expected <- pmax(c(32, 18) - alpha * spectral2, 0) /
    (c(32, 18) + beta * spectral2)
  for (fun in list(spca, rspca, robspca)) {
    fit <- fun(X, 2, alpha = alpha, beta = beta, penalty = "elasticnet",
               max_iter = 1000, tol = 1e-13, verbose = FALSE)
    expect_equal(abs(fit$loadings), rbind(diag(expected), c(0, 0)),
                 tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(fit$alpha, alpha * spectral2, tolerance = 1e-12)
    expect_equal(fit$beta, beta * spectral2, tolerance = 1e-12)
    expect_identical(fit$penalty, "elasticnet")
  }
})

test_that("l0 hard thresholding solves the one-component orthogonal problem", {
  X <- rbind(diag(c(4, 3, 2)), -diag(c(4, 3, 2)))
  # The active optimum is b = 1/(1 + beta); it is preferable to zero
  # precisely when alpha < 1/(2 * (1 + beta)).
  for (fun in list(spca, rspca, robspca)) {
    active <- fun(X, 1, alpha = 0.1, beta = 0.5, penalty = "l0",
                  tol = 1e-12, verbose = FALSE)
    inactive <- fun(X, 1, alpha = 0.5, beta = 0.5, penalty = "l0",
                    tol = 1e-12, verbose = FALSE)
    expect_equal(abs(active$loadings[, 1]), c(2/3, 0, 0), tolerance = 1e-10,
                 ignore_attr = TRUE)
    expect_equal(inactive$loadings, matrix(0, 3, 1), ignore_attr = TRUE)
    expect_identical(active$penalty, "l0")
    expect_true(all(diff(active$objective) <= 1e-10))
    expect_true(all(diff(inactive$objective) <= 1e-10))
    expected <- 0.5 * sum((X - active$scores %*% t(active$transform))^2) +
      3.2 * sum(active$loadings != 0) + 8 * sum(active$loadings^2)
    expect_equal(tail(active$objective, 1), expected, tolerance = 1e-10)
  }
})

test_that("the default penalty remains elastic net", {
  set.seed(313)
  X <- matrix(rnorm(20 * 6), 20, 6)
  for (fun in list(spca, rspca, robspca)) {
    set.seed(410)
    implicit <- fun(X, 2, max_iter = 20, verbose = FALSE)
    set.seed(410)
    explicit <- fun(X, 2, penalty = "elasticnet", max_iter = 20, verbose = FALSE)
    expect_equal(implicit$loadings, explicit$loadings)
    expect_equal(implicit$objective, explicit$objective)
  }
})

test_that("a loading warm start retains the requested mathematical convention", {
  set.seed(640)
  X <- matrix(rnorm(35 * 7), 35, 7)
  original <- spca(X, 2, alpha = 0.02, beta = 0.05,
                   max_iter = 1000, tol = 1e-12, verbose = FALSE)
  continued <- spca(X, 2, alpha = 0.02, beta = 0.05, B_init = original$loadings,
                    max_iter = 20, tol = 1e-12, verbose = FALSE)
  expect_lte(tail(continued$objective, 1), tail(original$objective, 1) * (1 + 1e-10))
  expect_equal(continued$scores,
               sweep(X, 2, continued$center, "-") %*% continued$loadings,
               tolerance = 1e-12, ignore_attr = TRUE)

  reference <- prcomp(X, rank. = 2)
  signed <- spca(X, 2, alpha = 0, beta = 0, B_init = -reference$rotation,
                verbose = FALSE)
  expect_equal(signed$scores %*% t(signed$transform),
               reference$x %*% t(reference$rotation), tolerance = 1e-9,
               ignore_attr = TRUE)
})

test_that("penalty and warm-start validation reject ambiguous inputs", {
  X <- matrix(seq_len(30), 6, 5)
  for (fun in list(spca, rspca, robspca)) {
    expect_error(fun(X, 2, penalty = "other", verbose = FALSE))
    expect_error(fun(X, 2, penalty = NA_character_, verbose = FALSE))
    for (initial in list(matrix(0, 4, 2), matrix(0, 5, 3), rep(0, 10),
                         matrix(NA_real_, 5, 2), matrix(Inf, 5, 2),
                         matrix(1i, 5, 2))) {
      expect_error(fun(X, 2, B_init = initial, verbose = FALSE))
    }
  }
})
