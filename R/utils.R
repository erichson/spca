# Input and numerical helpers shared by the three public solvers.
.spca_scalar <- function(x, name, lower = 0, integer = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || x < lower || (integer && x != floor(x))) {
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " >= ", lower, ".", call. = FALSE)
  }
  x
}

.spca_flag <- function(x, name) {
  if (length(x) != 1L || is.na(x) ||
      !(is.logical(x) || is.numeric(x)) || !(x %in% c(0, 1))) {
    stop(name, " must be TRUE or FALSE (or 1 or 0).", call. = FALSE)
  }
  as.logical(x)
}

.spca_matrix <- function(X, omit = TRUE, sparse = FALSE, min_rows = 2L) {
  if (inherits(X, "sparseMatrix") && sparse) {
    X <- methods::as(methods::as(methods::as(X, "generalMatrix"),
                               "dMatrix"), "CsparseMatrix")
    missing <- is.na(X@x)
    if (any(missing)) {
      if (!omit) stop("newdata must not contain missing values.", call. = FALSE)
      warning("Missing values are omitted: na.omit(X).", call. = FALSE)
      rows <- unique(X@i[missing]) + 1L
      X <- X[-rows, , drop = FALSE]
    }
    finite <- all(is.finite(X@x))
  } else {
    X <- as.matrix(X)
    if (length(dim(X)) != 2L || is.complex(X) ||
        !(is.numeric(X) || is.logical(X))) {
      stop("X must be a real numeric matrix or data frame.", call. = FALSE)
    }
    if (anyNA(X)) {
      if (!omit) stop("newdata must not contain missing values.", call. = FALSE)
      warning("Missing values are omitted: na.omit(X).", call. = FALSE)
      X <- X[rowSums(is.na(X)) == 0L, , drop = FALSE]
    }
    storage.mode(X) <- "double"
    finite <- all(is.finite(X))
  }
  if (nrow(X) < min_rows || ncol(X) < 1L) {
    stop("X must have at least ", min_rows,
         " rows and one column after omitting missing rows.", call. = FALSE)
  }
  if (!finite) stop("X must contain only finite values.", call. = FALSE)
  X
}

.spca_controls <- function(k, n, p, alpha, beta, max_iter, tol, verbose) {
  if (is.null(k)) k <- min(n, p)
  k <- min(.spca_scalar(k, "k", 1, TRUE), n, p)
  if (.spca_scalar(max_iter, "max_iter", 1, TRUE) > .Machine$integer.max) {
    stop("max_iter is too large.", call. = FALSE)
  }
  list(k = as.integer(k), alpha = .spca_scalar(alpha, "alpha"),
       beta = .spca_scalar(beta, "beta"),
       max_iter = as.integer(.spca_scalar(max_iter, "max_iter", 1, TRUE)),
       tol = .spca_scalar(tol, "tol"), verbose = .spca_flag(verbose, "verbose"))
}

.spca_start <- function(B_init, p, k, initial) {
  if (is.null(B_init)) return(initial)
  B <- as.matrix(B_init)
  if (!is.numeric(B) || is.complex(B) || !identical(dim(B), c(p, k)) ||
      any(!is.finite(B))) {
    stop("B_init must be a finite real matrix with p rows and k columns.",
         call. = FALSE)
  }
  unname(B)
}

.spca_prox <- function(z, threshold, penalty) {
  if (penalty == "l0") {
    z[abs(z) <= sqrt(2 * threshold)] <- 0
    z
  } else {
    sign(z) * pmax(abs(z) - threshold, 0)
  }
}

.spca_penalty <- function(B, alpha, beta, penalty) {
  alpha * if (penalty == "l0") sum(B != 0) else sum(abs(B))
}

.spca_objective <- function(loss, B, alpha, beta, penalty) {
  loss + .spca_penalty(B, alpha, beta, penalty) + 0.5 * beta * sum(B^2)
}

.spca_progress <- function(previous, current, tol, evaluation_scale = current) {
  if (!is.finite(current)) {
    stop("The objective is not finite; rescale X or reduce the penalties.",
         call. = FALSE)
  }
  if (is.null(previous)) return(list(improvement = Inf, converged = FALSE,
                                    stop = FALSE))
  delta <- previous - current
  roundoff <- 64 * .Machine$double.eps *
    max(abs(previous), abs(current), abs(evaluation_scale))
  if (delta < -roundoff) {
    warning("The objective increased; stopping without convergence.", call. = FALSE)
    return(list(improvement = delta / max(abs(current), .Machine$double.xmin),
                converged = FALSE, stop = TRUE))
  }
  improvement <- max(0, delta) / max(abs(current), .Machine$double.xmin)
  list(improvement = improvement, converged = improvement <= tol,
       stop = improvement <= tol)
}

.spca_trace <- function(i, objective, improvement, verbose) {
  if (verbose && (i == 1L || i %% 10L == 0L)) {
    cat(sprintf("Iteration: %4d, Objective: %1.5e, Relative improvement %1.5e\n",
                i, objective, improvement))
  }
}

.spca_result <- function(fit, prepared, class, alpha, beta, penalty) {
  scores <- prepared$forward(fit$B)
  rownames(scores) <- prepared$row_names
  # dqrdc2 only moves numerically dependent columns to the end. Map the
  # retained diagonal energies back to the original component order.
  qr_scores <- qr(scores, LAPACK = FALSE)
  adjusted <- numeric(ncol(scores))
  if (qr_scores$rank > 0L) {
    keep <- seq_len(qr_scores$rank)
    adjusted[qr_scores$pivot[keep]] <- diag(qr.R(qr_scores))[keep]^2
  }
  eigenvalues <- adjusted / (prepared$n - 1)
  out <- list(loadings = fit$B, transform = fit$A, scores = scores,
              eigenvalues = eigenvalues, center = prepared$center,
              scale = prepared$scale)
  if (!is.null(fit$S)) out$sparse <- fit$S
  out$objective <- fit$objective
  out$sdev <- sqrt(eigenvalues)
  out$var <- prepared$energy / (prepared$n - 1)
  out$score_variance <- colSums(scores^2) / (prepared$n - 1)
  out$iterations <- length(fit$objective)
  out$converged <- fit$converged
  out$improvement <- fit$improvement
  out$alpha <- alpha
  out$beta <- beta
  out$penalty <- penalty
  out$feature_names <- prepared$feature_names
  class(out) <- class
  out
}
