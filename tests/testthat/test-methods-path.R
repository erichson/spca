test_that("prediction reuses training preprocessing and loading coefficients", {
  X <- as.matrix(iris[, 1:4])
  for (fun in list(spca, rspca, robspca)) {
    for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
      fit <- fun(X[1:100, ], 2, center = center, scale = scale, verbose = FALSE)
      explicit <- X
      if (center) explicit <- sweep(explicit, 2, fit$center, "-")
      if (scale) explicit <- sweep(explicit, 2, fit$scale, "/")
      expected <- explicit %*% fit$loadings
      expect_equal(predict(fit), fit$scores)
      expect_equal(predict(fit, X), expected, tolerance = 1e-12)
      expect_equal(predict(fit, X[, 4:1]), expected, tolerance = 1e-12)
      expect_equal(predict(fit, Matrix::Matrix(X, sparse = TRUE)), expected,
                   tolerance = 1e-12)
      expect_equal(predict(fit, X[1, , drop = FALSE]), expected[1, , drop = FALSE],
                   tolerance = 1e-12)
      expect_identical(dim(predict(fit, X[FALSE, , drop = FALSE])), c(0L, 2L))
      expect_identical(dim(summary(fit)), c(4L, 2L))
    }
  }
})

test_that("prediction rejects mismatched or missing feature values", {
  X <- as.matrix(iris[, 1:4])
  fit <- spca(X, 2, verbose = FALSE)
  expect_error(predict(fit, X[, 1:3]), "columns")
  wrong <- X
  colnames(wrong)[1] <- "not-a-feature"
  expect_error(predict(fit, wrong), "column names")
  colnames(wrong) <- rep("feature", 4)
  expect_error(predict(fit, wrong), "column names")
  wrong <- X
  wrong[1, 1] <- NA
  expect_error(predict(fit, wrong), "missing")
  expect_error(predict(fit, Matrix::Matrix(wrong, sparse = TRUE)), "missing")
})

test_that("a path shares a factorization and follows the given alpha order", {
  set.seed(74)
  X <- matrix(rnorm(60 * 8), 60, 8)
  alpha <- c(0.003, 0.01, 0.004)
  path <- spca_path(X, 2, alpha = alpha, max_iter = 60, tol = 1e-9,
                    verbose = FALSE)
  expect_s3_class(path, "spca_path")
  expect_identical(path$alpha, alpha)
  first <- spca(X, 2, alpha = alpha[1], max_iter = 60, tol = 1e-9,
                verbose = FALSE)
  second <- spca(X, 2, alpha = alpha[2], B_init = first$loadings,
                 max_iter = 60, tol = 1e-9, verbose = FALSE)
  expect_equal(path$fits[[1]]$loadings, first$loadings, tolerance = 1e-12)
  expect_equal(path$fits[[2]]$loadings, second$loadings, tolerance = 1e-12)
  expect_equal(unname(path$nonzero[2, ]), unname(colSums(second$loadings != 0)))
  expect_equal(path$selected_features[2], sum(rowSums(second$loadings != 0) > 0))
  expect_identical(path$iterations, vapply(path$fits, `[[`, integer(1), "iterations"))
  expect_output(print(path), "selected_features")
  expect_error(spca_path(X, 2, alpha = numeric()), "alpha")
  expect_error(spca_path(X, 2, alpha = c(0.1, NA)), "alpha")
})

test_that("randomized paths use one sketch and respect the R RNG", {
  set.seed(9)
  X <- matrix(rnorm(80 * 25), 80, 25)
  set.seed(910)
  before <- .Random.seed
  first <- rspca(X, 2, o = 3, q = 1, max_iter = 2, verbose = FALSE)
  after_one <- .Random.seed
  expect_false(identical(before, after_one))
  set.seed(910)
  path <- spca_path(X, 2, method = "rspca", o = 3, q = 1,
                    alpha = c(1e-4, 1e-3), max_iter = 2, verbose = FALSE)
  expect_identical(.Random.seed, after_one)
  expect_equal(path$fits[[1]]$loadings, first$loadings, tolerance = 1e-12)
  seed <- .Random.seed
  spca(X, 2, max_iter = 2, verbose = FALSE)
  robspca(X, 2, max_iter = 2, verbose = FALSE)
  expect_identical(.Random.seed, seed)
})
