test_that("implicit centered and scaled products match explicit dense products", {
  set.seed(333)
  base <- matrix(rnorm(24 * 9), 24, 9)
  base[abs(base) < 0.8] <- 0
  base[, 9] <- 0
  matrices <- list(Matrix::Matrix(base, sparse = TRUE),
                   methods::as(Matrix::Matrix(base, sparse = TRUE), "TsparseMatrix"),
                   methods::as(Matrix::Matrix(base, sparse = TRUE), "RsparseMatrix"),
                   Matrix::Diagonal(x = seq_len(9)),
                   Matrix::forceSymmetric(Matrix::Matrix(base[1:9, ], sparse = TRUE)),
                   Matrix::triu(Matrix::Matrix(base[1:9, ], sparse = TRUE)),
                   Matrix::Matrix(base != 0, sparse = TRUE),
                   methods::as(Matrix::Matrix(base != 0, sparse = TRUE), "nMatrix"))
  for (X in matrices) {
    for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
      p <- sparsepca:::.spca_preprocess(X, center, scale, sparse = TRUE)
      dense <- sparsepca:::.spca_preprocess(as.matrix(X), center, scale)
      right <- matrix(rnorm(ncol(X) * 3), ncol(X), 3)
      left <- matrix(rnorm(nrow(X) * 3), nrow(X), 3)
      expect_equal(p$forward(right), dense$forward(right), tolerance = 1e-12,
                   ignore_attr = TRUE)
      expect_equal(p$adjoint(left), dense$adjoint(left), tolerance = 1e-12,
                   ignore_attr = TRUE)
      expect_equal(p$energy, dense$energy, tolerance = 1e-12)
      expect_equal(sum(p$forward(right) * left), sum(right * p$adjoint(left)),
                   tolerance = 1e-12)
    }
  }
})

test_that("fully stored sparse columns retain accuracy under large offsets", {
  set.seed(54)
  X <- matrix(1e12 + rnorm(50 * 7), 50, 7)
  sparse <- Matrix::Matrix(X, sparse = TRUE)
  for (scale in c(FALSE, TRUE)) {
    sparse_p <- sparsepca:::.spca_preprocess(sparse, TRUE, scale, sparse = TRUE)
    dense_p <- sparsepca:::.spca_preprocess(X, TRUE, scale)
    W <- matrix(rnorm(7 * 3), 7, 3)
    Y <- matrix(rnorm(50 * 3), 50, 3)
    # Matrix versions differ in their representation of empty dimnames.
    # Retain checks of matrix shape and values while ignoring only names.
    expect_equal(unname(sparse_p$forward(W)), unname(dense_p$forward(W)),
                 tolerance = 1e-12)
    expect_equal(unname(sparse_p$adjoint(Y)), unname(dense_p$adjoint(Y)),
                 tolerance = 1e-12)
    expect_equal(sparse_p$energy, dense_p$energy, tolerance = 1e-12)
  }
})

test_that("sparse and dense fits have the same statistical meaning", {
  set.seed(417)
  X <- matrix(rnorm(40 * 12), 40, 12)
  X[abs(X) < 1] <- 0
  sparse <- Matrix::Matrix(X, sparse = TRUE)
  for (fun in list(spca, rspca, robspca)) {
    for (scale in c(FALSE, TRUE)) {
      set.seed(55)
      dense_fit <- fun(X, 3, scale = scale, alpha = 0.005,
                       max_iter = 500, tol = 1e-9, verbose = FALSE)
      set.seed(55)
      sparse_fit <- fun(sparse, 3, scale = scale, alpha = 0.005,
                        max_iter = 500, tol = 1e-9, verbose = FALSE)
      expect_equal(sparse_fit$scores %*% t(sparse_fit$transform),
                   dense_fit$scores %*% t(dense_fit$transform), tolerance = 1e-7)
      expect_equal(sparse_fit$var, dense_fit$var, tolerance = 1e-12)
      expect_equal(sparse_fit$eigenvalues, dense_fit$eigenvalues, tolerance = 1e-7)
    }
  }
})

test_that("randomized sparse fits avoid an n by p dense allocation", {
  skip_if_not(capabilities("profmem"))
  set.seed(99)
  X <- Matrix::rsparsematrix(700, 400, density = 0.01)
  path <- tempfile()
  on.exit(unlink(path))
  Rprofmem(path, threshold = 8 * prod(dim(X)) - 1000)
  tryCatch(rspca(X, 3, o = 4, max_iter = 2, verbose = FALSE),
           finally = Rprofmem(NULL))
  allocations <- readLines(path)
  sizes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
  expect_false(any(sizes >= 8 * prod(dim(X)), na.rm = TRUE))
})

test_that("sparse missing rows and structural unit diagonals are respected", {
  X <- Matrix::sparseMatrix(i = c(1, 3, 4), j = c(1, 2, 3),
                            x = c(2, NA, 4), dims = c(6, 4))
  expect_warning(fit <- rspca(X, 2, verbose = FALSE), "Missing")
  expect_equal(nrow(fit$scores), 5)
  unit <- methods::as(Matrix::Diagonal(5), "triangularMatrix")
  set.seed(6)
  a <- rspca(unit, 2, alpha = 0, beta = 0, verbose = FALSE)
  set.seed(6)
  b <- rspca(as.matrix(unit), 2, alpha = 0, beta = 0, verbose = FALSE)
  # Four tied directions permit different rank-two bases. Check the optimal
  # residual and energies, not equality of a particular arbitrary basis.
  centered <- scale(as.matrix(unit), scale = FALSE)
  for (fit in list(a, b)) {
    expect_equal(sum((centered - tcrossprod(fit$scores, fit$transform))^2), 2,
                 tolerance = 1e-10)
    expect_equal(fit$eigenvalues, c(0.25, 0.25), tolerance = 1e-10)
  }
})
