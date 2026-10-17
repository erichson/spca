test_that("a genuinely compressed sparse sketch recovers a sparse low-rank matrix", {
  # Disjoint feature blocks and pairs of opposite observations give an
  # exactly rank-four matrix with 1.67% nonzeros. The seven sketch rows are
  # far fewer than min(n, p) = 72, so this cannot pass via a full-width sketch.
  n <- 120L
  p <- 72L
  k <- 4L
  left <- matrix(0, n, k)
  right <- matrix(0, p, k)
  for (j in seq_len(k)) {
    left[(2L * j - 1L):(2L * j), j] <- c(1, -1) / sqrt(2)
    rows <- seq.int((j - 1L) * 18L + 1L, j * 18L)
    right[rows, j] <- seq_along(rows) / sqrt(sum(seq_along(rows)^2))
  }
  X <- left %*% diag(c(30, 15, 7, 3)) %*% t(right)
  sparse <- Matrix::Matrix(X, sparse = TRUE)
  for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
    # Base R supplies the explicit reference, independently of either
    # sparsepca preprocessing path. Every column has a nonzero scale.
    target <- base::scale(X, center = center, scale = scale)
    for (q in c(0, 2)) {
      set.seed(21)
      fit <- rspca(sparse, k, alpha = 0, beta = 0, center = center,
                    scale = scale, o = 3, q = q, verbose = FALSE)
      reconstruction <- tcrossprod(fit$scores, fit$transform)
      expect_lt(norm(target - reconstruction, "F") / norm(target, "F"), 1e-10)
      expect_equal(crossprod(fit$transform), diag(k), tolerance = 1e-10)
      expect_equal(fit$var, sum(target^2) / (n - 1), tolerance = 1e-12)
    }
  }
})

test_that("power iterations improve a compressed sparse fit toward the spectral tail", {
  n <- 120L
  p <- 72L
  k <- 4L
  spectrum <- exp(seq(log(30), log(0.03), length.out = p))
  sparse <- Matrix::sparseMatrix(i = seq_len(p), j = seq_len(p),
                                  x = spectrum, dims = c(n, p))
  for (center in c(FALSE, TRUE)) {
    # This full-rank input cannot be reconstructed exactly at rank four.
    # The SVD tail gives the optimal attainable Frobenius error, including
    # the effect of centering when that changes the known diagonal spectrum.
    target <- base::scale(as.matrix(sparse), center = center, scale = FALSE)
    singular <- svd(target, nu = 0, nv = 0)$d
    optimal_error <- sum(singular[-seq_len(k)]^2)
    for (seed in c(8, 21, 235)) {
      errors <- vapply(c(0, 4), function(q) {
        set.seed(seed)
        fit <- rspca(sparse, k, alpha = 0, beta = 0, center = center,
                      o = 3, q = q, verbose = FALSE)
        sum((target - tcrossprod(fit$scores, fit$transform))^2)
      }, numeric(1))
      # Observed ratios were <= 1.011 across these seeds; 1.03 leaves a
      # margin for numerical backends. No basis or sign is prescribed.
      expect_lte(errors[2], 1.03 * optimal_error)
      expect_gte(errors[2], (1 - 1e-10) * optimal_error)
      expect_lt(errors[2], 0.99 * errors[1])
    }
  }
})

test_that("mixed sparse offsets and scales match an explicit independent operator", {
  n <- 40L
  X <- cbind(full = 1e12 + rep(c(-3, -1, 1, 3), length.out = n),
             incomplete = c(0, 1e12 + seq_len(n - 1L)),
             sparse = rep(c(0, -2, 0, 3), length.out = n),
             tiny = rep(c(0, -1, 0, 1), length.out = n) * 1e-10,
             zero = rep(0, n))
  sparse <- Matrix::Matrix(X, sparse = TRUE)
  original <- sparse
  Y <- matrix(sin(seq_len(n * 3L)), n, 3L)
  for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
    target <- base::scale(X, center = center, scale = FALSE)
    scales <- vapply(seq_len(ncol(X)), function(j) {
      sqrt(sum(target[, j]^2) / (n - 1L))
    }, numeric(1))
    # Preserve the documented small-scale convention, including zero columns.
    scales[scales < 1e-8] <- 1
    if (scale) target <- sweep(target, 2, scales, "/")
    prepared <- sparsepca:::.spca_preprocess(sparse, center, scale, sparse = TRUE)
    forward <- prepared$forward(diag(ncol(X)))
    adjoint <- prepared$adjoint(Y)
    reference_adjoint <- crossprod(target, Y)
    # Compare columns separately: the large incomplete column must not hide
    # cancellation errors in the fully stored, explicitly centered column.
    for (j in seq_len(ncol(X))) {
      expect_equal(unname(forward[, j]), unname(target[, j]), tolerance = 1e-12)
      expect_equal(unname(adjoint[j, ]), unname(reference_adjoint[j, ]),
                   tolerance = 1e-12)
    }
    expect_equal(prepared$energy, sum(target^2), tolerance = 1e-12)
  }
  expect_identical(sparse, original)
})
