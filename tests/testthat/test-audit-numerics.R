test_that("underflow cannot turn nonzero data into a converged zero fit", {
  # The unpenalized, one-variable problem has |B| = 1 at every nonzero
  # scale. These entries are representable, but their squares are not.
  # Returning B = 0 with objective 0 would falsely report an exact fit.
  X <- matrix(c(-1, 1), ncol = 1) * 1e-170
  for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
    for (fun in list(spca, rspca, robspca)) {
      expect_error(fun(X, 1, alpha = 0, beta = 0, center = center,
                       scale = scale, verbose = FALSE),
                   "squared data norm underflows to zero; rescale X")
    }
    expect_error(rspca(Matrix::Matrix(X, sparse = TRUE), 1,
                       alpha = 0, beta = 0, center = center, scale = scale,
                       verbose = FALSE),
                 "squared data norm underflows to zero; rescale X")
  }
})

test_that("zero-energy validation preserves genuine centered constants", {
  # Large constants become exactly zero after centering. The guard must
  # examine processed entries, not the original nonzero values.
  X <- matrix(1e170, 5, 2)
  for (input in list(X, Matrix::Matrix(X, sparse = TRUE))) {
    fit <- rspca(input, 1, scale = TRUE, verbose = FALSE)
    expect_true(fit$converged)
    expect_equal(fit$loadings, matrix(0, 2, 1))
    # Older Matrix versions retain empty dimnames on sparse products.
    expect_equal(unname(fit$scores), matrix(0, 5, 1))
    expect_identical(fit$objective, 0)
  }
})

test_that("penalty searches reject underflow before selecting a zero fit", {
  X <- matrix(c(-1, 1), ncol = 1) * 1e-170
  for (method in c("spca", "rspca")) {
    input <- if (method == "rspca") Matrix::Matrix(X, sparse = TRUE) else X
    expect_error(spca_path(input, 1, alpha = c(0, .01), method = method,
                           verbose = FALSE), "underflows to zero; rescale X")
    expect_error(spca_target(input, 1, n_features = 1, method = method,
                             verbose = FALSE), "underflows to zero; rescale X")
  }
})
