.spca_preprocess <- function(X, center, scale, sparse = FALSE) {
  center <- .spca_flag(center, "center")
  scale <- .spca_flag(scale, "scale")
  X <- .spca_matrix(X, sparse = sparse)
  n <- nrow(X)
  p <- ncol(X)
  row_names <- rownames(X)
  feature_names <- colnames(X)
  if (inherits(X, "sparseMatrix")) {
    counts <- diff(X@p)
    means <- Matrix::colMeans(X)
    full <- counts == n
    for (j in which(full)) {
      ids <- seq.int(X@p[j] + 1L, X@p[j + 1L])
      # Match dense base colMeans, including its extended-precision sum.
      means[j] <- colMeans(matrix(X@x[ids], ncol = 1L))
    }
    for (j in which(!full & !is.finite(means))) {
      ids <- seq.int(X@p[j] + 1L, X@p[j + 1L])
      means[j] <- mean(X@x[ids]) * (counts[j] / n)
    }
    offset <- if (center) means else rep(0, p)
    # Accumulate stored deviations and implicit zeros separately, avoiding
    # sum(x^2) - n*mean(x)^2 for almost constant columns.
    deviations <- X
    deviations@x <- (X@x - rep(offset, counts))^2
    squared <- Matrix::colSums(deviations)
    incomplete <- counts < n
    squared[incomplete] <- squared[incomplete] +
      (n - counts[incomplete]) * offset[incomplete]^2
    scales <- if (scale) sqrt(squared / (n - 1)) else rep(1, p)
    scales[scales < 1e-8] <- 1
    if (any(!is.finite(scales))) {
      stop("The column scales are not finite; rescale X.", call. = FALSE)
    }
    energy <- sum(squared / scales^2)
    stored_center <- if (center) means else FALSE
    stored_scale <- if (scale) scales else FALSE
    products <- .spca_sparse_products(X, stored_center, stored_scale)
    data <- NULL
  } else {
    stored_center <- if (center) colMeans(X) else FALSE
    if (center) X <- sweep(X, 2L, stored_center, "-")
    stored_scale <- if (scale) sqrt(colSums(X^2) / (n - 1)) else FALSE
    if (scale) {
      stored_scale[stored_scale < 1e-8] <- 1
      if (any(!is.finite(stored_scale))) {
        stop("The column scales are not finite; rescale X.", call. = FALSE)
      }
      X <- sweep(X, 2L, stored_scale, "/")
    }
    energy <- sum(X^2)
    products <- list(forward = function(B) X %*% B,
                     adjoint = function(Y) crossprod(X, Y))
    data <- X
  }
  if (!is.finite(energy)) {
    stop("The squared data norm is not finite; rescale X.", call. = FALSE)
  }
  list(data = data, n = n, p = p, energy = energy,
       center = stored_center, scale = stored_scale,
       forward = products$forward, adjoint = products$adjoint,
       row_names = row_names, feature_names = feature_names)
}

.spca_sparse_products <- function(X, center, scale) {
  # X is a canonical real general Csparse matrix, including any formerly
  # implicit diagonal/symmetric entries. Fully stored columns can be centered
  # in place; this protects products with large offsets against cancellation.
  counts <- diff(X@p)
  shift <- if (identical(center, FALSE)) rep(0, ncol(X)) else center
  full <- counts == nrow(X)
  explicit <- shift
  explicit[!full] <- 0
  X@x <- X@x - rep(explicit, counts)
  shift[full] <- 0
  if (!identical(scale, FALSE)) {
    X@x <- X@x / rep(scale, counts)
    shift <- shift / scale
  }
  X@factors <- list()
  centered <- any(shift != 0)
  list(
    forward = function(B) {
      Y <- as.matrix(X %*% B)
      if (centered) Y <- sweep(Y, 2L, as.vector(crossprod(shift, B)), "-")
      Y
    },
    adjoint = function(Y) {
      B <- as.matrix(Matrix::crossprod(X, Y))
      if (centered) B <- B - tcrossprod(shift, colSums(Y))
      B
    }
  )
}

.spca_sparse_sketch <- function(prepared, l, q) {
  omega <- matrix(stats::rnorm(prepared$p * l), prepared$p, l)
  Y <- prepared$forward(omega)
  if (q > 0L) {
    for (i in seq_len(q)) {
      Q <- qr.Q(qr(Y), complete = FALSE)
      Z <- qr.Q(qr(prepared$adjoint(Q)), complete = FALSE)
      Y <- prepared$forward(Z)
    }
  }
  Q <- qr.Q(qr(Y), complete = FALSE)
  t(prepared$adjoint(Q))
}

.spca_factor <- function(prepared, k, randomized = FALSE, o = 20, q = 2) {
  if (randomized) {
    o <- .spca_scalar(o, "o", 0, TRUE)
    q <- .spca_scalar(q, "q", 0, TRUE)
  }
  if (prepared$energy == 0) {
    return(list(C = matrix(0, 1L, prepared$p), norm2 = 0,
                initial = diag(1, prepared$p, k)))
  }
  X <- prepared$data
  if (randomized) {
    l <- min(k + o, prepared$n, prepared$p)
    X <- if (is.null(X)) .spca_sparse_sketch(prepared, l, q) else
      rsvd::rqb(X, k = l, p = 0, q = q)$B
  }
  decomp <- svd(X, nu = 0, nv = min(dim(X)))
  norm2 <- decomp$d[1L]^2
  list(C = decomp$d * t(decomp$v), norm2 = norm2,
       initial = decomp$v[, seq_len(k), drop = FALSE])
}
