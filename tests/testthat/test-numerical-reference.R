test_that("unpenalized SPCA matches PCA with the same preprocessing", {
  set.seed(119)
  for (dims in list(c(31, 8), c(8, 31))) {
    X <- matrix(rnorm(prod(dims)), dims[1], dims[2])
    X <- sweep(X, 2, seq_len(ncol(X))/3, "+")
    for (center in c(TRUE, FALSE)) {
      for (scale in c(TRUE, FALSE)) {
        fit <- spca(X, 3, alpha = 0, beta = 0, center = center, scale = scale,
                    verbose = FALSE)
        reference <- prcomp(X, rank. = 3, center = center, scale. = scale)
        processed <- scale(X, center = center, scale = scale)
        reconstructed <- fit$scores %*% t(fit$transform)
        optimal <- reference$x %*% t(reference$rotation)
        expect_equal(reconstructed, optimal, tolerance = 1e-9, ignore_attr = TRUE)
        expect_equal(fit$sdev, reference$sdev[1:3], tolerance = 1e-9,
                     ignore_attr = TRUE)
        expect_equal(crossprod(fit$transform), diag(3), tolerance = 1e-12,
                     ignore_attr = TRUE)
        expect_equal(fit$scores, processed %*% fit$loadings, tolerance = 1e-12,
                     ignore_attr = TRUE)
        expect_equal(sum((processed - reconstructed)^2),
                     sum((processed - optimal)^2), tolerance = 1e-9)
        expect_equal(fit$var, sum(processed^2)/(nrow(X) - 1), tolerance = 1e-12)
      }
    }
  }
})

test_that("randomized PCA attains the exact rank and a controlled spectral tail", {
  set.seed(492)
  left <- qr.Q(qr(matrix(rnorm(80 * 50), 80, 50)))
  right <- qr.Q(qr(matrix(rnorm(50 * 50), 50, 50)))
  exact <- sweep(left[, 1:5, drop = FALSE], 2, c(20, 15, 10, 5, 2), "*") %*%
    t(right[, 1:5, drop = FALSE])
  for (seed in c(8, 38)) {
    set.seed(seed)
    fit <- rspca(exact, 5, alpha = 0, beta = 0, center = FALSE,
                 o = 5, q = 1, verbose = FALSE)
    relative <- norm(exact - fit$scores %*% t(fit$transform), "F") / norm(exact, "F")
    expect_lt(relative, 1e-10)
    expect_equal(crossprod(fit$transform), diag(5), tolerance = 1e-12,
                 ignore_attr = TRUE)
  }

  spectrum <- 30 * seq_len(50)^(-0.85)
  X <- sweep(left, 2, spectrum, "*") %*% t(right)
  optimal_error2 <- sum(spectrum[-seq_len(6)]^2)
  for (seed in c(19, 40)) {
    set.seed(seed)
    fit <- rspca(X, 6, alpha = 0, beta = 0, center = FALSE,
                 o = 8, q = 3, verbose = FALSE)
    observed_error2 <- sum((X - fit$scores %*% t(fit$transform))^2)
    expect_gte(observed_error2, optimal_error2 * (1 - 1e-10))
    expect_lt(observed_error2/optimal_error2, 1.08)
  }
})

test_that("repeated singular values are checked as a subspace", {
  set.seed(499)
  left <- qr.Q(qr(matrix(rnorm(25 * 5), 25, 5)))
  right <- qr.Q(qr(matrix(rnorm(9 * 5), 9, 5)))
  X <- sweep(left, 2, c(10, 10, 10, 3, 1), "*") %*% t(right)
  fit <- spca(X, 3, alpha = 0, beta = 0, center = FALSE, verbose = FALSE)
  expect_equal(tcrossprod(fit$loadings),
               tcrossprod(right[, 1:3, drop = FALSE]),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(sum((X - fit$scores %*% t(fit$transform))^2), 10,
               tolerance = 1e-10)
})

test_that("rank-deficient data preserve requested dimensions and finite diagnostics", {
  X <- outer(seq_len(20) - 10.5, c(3, 0, -1, 2, 0, 5, 1))
  for (fun in list(spca, rspca, robspca)) {
    fit <- fun(X, 3, alpha = 0, beta = 0, verbose = FALSE)
    expect_identical(dim(fit$loadings), c(7L, 3L))
    expect_true(all(is.finite(unlist(fit[c("eigenvalues", "sdev", "objective")]))))
    expect_lt(norm(X - fit$scores %*% t(fit$transform), "F")/norm(X, "F"), 1e-10)
    expect_lt(sum(fit$eigenvalues[-1]), 1e-10 * fit$eigenvalues[1])
  }
})

test_that("variance diagnostics use returned scores and account for correlation", {
  # An independent sequential regression reference expresses the adjusted
  # variance without reusing the package's QR implementation.
  adjusted_energy <- function(scores) {
    vapply(seq_len(ncol(scores)), function(j) {
      residual <- scores[, j]
      if (j > 1L) {
        previous <- scores[, seq_len(j - 1L), drop = FALSE]
        residual <- qr.resid(qr(previous), residual)
      }
      sum(residual^2)/(nrow(scores) - 1)
    }, numeric(1))
  }
  set.seed(18)
  X <- matrix(rnorm(50 * 9), 50, 9)
  X[, 5:9] <- X[, 5:9] + 0.8 * X[, 1:5]
  for (center in c(TRUE, FALSE)) {
    fit <- spca(X + 2, 3, alpha = 0.03, beta = 0.08,
                center = center, max_iter = 500, tol = 1e-9, verbose = FALSE)
    expect_equal(fit$eigenvalues, adjusted_energy(fit$scores), tolerance = 1e-9,
                 ignore_attr = TRUE)
    expect_equal(fit$sdev^2, fit$eigenvalues, tolerance = 1e-12,
                 ignore_attr = TRUE)
    expect_equal(fit$score_variance, colSums(fit$scores^2)/(nrow(X) - 1),
                 tolerance = 1e-12, ignore_attr = TRUE)
    expect_lte(sum(fit$eigenvalues), sum(fit$score_variance) * (1 + 1e-10))
    processed <- scale(X + 2, center = center, scale = FALSE)
    expect_equal(fit$var, sum(processed^2)/(nrow(X) - 1), tolerance = 1e-12)
    expect_equal(unname(summary(fit)[3, ]), round(fit$eigenvalues/fit$var, 3),
                 ignore_attr = TRUE)
  }
})

test_that("ridge shrinkage is retained and its score variance is measured correctly", {
  X <- rbind(diag(c(4, 3, 2)), -diag(c(4, 3, 2)))
  for (fun in list(spca, rspca, robspca)) {
    fit <- fun(X, 1, alpha = 0, beta = 0.5, verbose = FALSE)
    expect_equal(abs(fit$loadings[, 1]), c(2/3, 0, 0), tolerance = 1e-10,
                 ignore_attr = TRUE)
    expect_equal(fit$eigenvalues, (32/5) * (2/3)^2, tolerance = 1e-10,
                 ignore_attr = TRUE)
  }
})
