# The right SVD factor C = D V' has the same Gram matrix as the data (or
# sketch). Iterations use C times thin factors; they do not form a p by p
# covariance matrix or a full residual on every step.
.spca_loss <- function(C, CA, CB, A, energy) {
  value <- energy - 2 * sum(CA * CB) + sum(CB^2)
  # Near an exact reconstruction the difference of large energies is unsafe.
  # Recompute the actual residual in bounded row blocks instead of clipping
  # a potentially inaccurate objective to zero.
  if (value <= 1e-8 * energy) {
    value <- 0
    block <- max(1L, floor(1e6 / ncol(C)))
    for (start in seq.int(1L, nrow(C), by = block)) {
      rows <- seq.int(start, min(nrow(C), start + block - 1L))
      residual <- C[rows, , drop = FALSE] -
        tcrossprod(CB[rows, , drop = FALSE], A)
      value <- value + sum(residual^2)
    }
  }
  0.5 * value
}

.spca_solve <- function(factor, controls, penalty, B_init = NULL) {
  C <- factor$C
  k <- controls$k
  B <- .spca_start(B_init, ncol(C), k, factor$initial)
  if (factor$norm2 == 0) {
    return(list(B = matrix(0, ncol(C), k), A = factor$initial,
                objective = 0, converged = TRUE, improvement = 0))
  }
  alpha <- controls$alpha * factor$norm2
  beta <- controls$beta * factor$norm2
  if (!is.finite(alpha) || !is.finite(beta)) {
    stop("The effective penalties are not finite; rescale X or the penalties.",
         call. = FALSE)
  }
  shrink <- 1 + controls$beta
  threshold <- controls$alpha / shrink
  objective <- numeric(controls$max_iter)
  energy <- sum(C^2)
  CB <- C %*% B
  progress <- list(improvement = Inf, converged = FALSE)
  for (i in seq_len(controls$max_iter)) {
    polar <- svd(crossprod(C, CB), nu = k, nv = k)
    A <- tcrossprod(polar$u, polar$v)
    CA <- C %*% A
    # Algebraically the legacy step 1/(||C||_2^2 + beta_eff), with the
    # ridge term in the smooth gradient; this form avoids a redundant copy.
    B <- .spca_prox((B + crossprod(C, CA - CB) / factor$norm2) / shrink,
                    threshold, penalty)
    CB <- C %*% B
    loss <- .spca_loss(C, CA, CB, A, energy)
    objective[i] <- .spca_objective(loss, B, alpha, beta, penalty)
    previous <- if (i > 1L) objective[i - 1L] else NULL
    progress <- .spca_progress(previous, objective[i], controls$tol,
                                evaluation_scale = energy)
    .spca_trace(i, objective[i], progress$improvement, controls$verbose)
    if (progress$stop) break
  }
  list(B = B, A = A, objective = objective[seq_len(i)],
       converged = progress$converged, improvement = progress$improvement)
}

.spca_robust_solve <- function(prepared, controls, gamma, penalty, B_init) {
  X <- prepared$data
  k <- controls$k
  init <- svd(X, nu = 0, nv = k)
  norm2 <- init$d[1L]^2
  B <- .spca_start(B_init, ncol(X), k, init$v[, seq_len(k), drop = FALSE])
  S <- matrix(0, nrow(X), ncol(X))
  if (norm2 == 0) {
    return(list(B = matrix(0, ncol(X), k), A = init$v,
                S = S, objective = 0, converged = TRUE, improvement = 0,
                norm2 = 0))
  }
  alpha <- controls$alpha * norm2
  beta <- controls$beta * norm2
  if (!is.finite(alpha) || !is.finite(beta)) {
    stop("The effective penalties are not finite; rescale X or the penalties.",
         call. = FALSE)
  }
  shrink <- 1 + controls$beta
  threshold <- controls$alpha / shrink
  objective <- numeric(controls$max_iter)
  XB <- X %*% B
  progress <- list(improvement = Inf, converged = FALSE)
  for (i in seq_len(controls$max_iter)) {
    XS <- X - S
    polar <- svd(crossprod(XS, XB), nu = k, nv = k)
    A <- tcrossprod(polar$u, polar$v)
    B <- .spca_prox((B + crossprod(X, XS %*% A - XB) / norm2) / shrink,
                    threshold, penalty)
    XB <- X %*% B
    R <- X - tcrossprod(XB, A)
    S <- sign(R) * pmax(abs(R) - gamma, 0)
    # R is the low-rank residual. The reconstruction loss must include S;
    # omitting it can stop after two iterations as soon as outliers appear.
    objective[i] <- .spca_objective(0.5 * sum((R - S)^2), B,
                                    alpha, beta, penalty) + gamma * sum(abs(S))
    previous <- if (i > 1L) objective[i - 1L] else NULL
    progress <- .spca_progress(previous, objective[i], controls$tol,
                                evaluation_scale = prepared$energy)
    .spca_trace(i, objective[i], progress$improvement, controls$verbose)
    if (progress$stop) break
  }
  list(B = B, A = A, S = S, objective = objective[seq_len(i)],
       converged = progress$converged, improvement = progress$improvement,
       norm2 = norm2)
}
