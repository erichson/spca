#' @rdname robspca
#' @export
robspca <- function(X, k = NULL, alpha = 1e-4, beta = 1e-4, gamma = 100, center = TRUE, scale = FALSE, max_iter = 1000, tol = 1e-5, verbose = TRUE, penalty = c("elasticnet", "l0"), B_init = NULL) {
  UseMethod("robspca")
}

#' @export
robspca.default <- function(X, k = NULL, alpha = 1e-4, beta = 1e-4, gamma = 100, center = TRUE, scale = FALSE, max_iter = 1000, tol = 1e-5, verbose = TRUE, penalty = c("elasticnet", "l0"), B_init = NULL) {
  penalty <- match.arg(penalty)
  prepared <- .spca_preprocess(X, center, scale, sparse = FALSE)
  controls <- .spca_controls(k, prepared$n, prepared$p, alpha, beta,
                              max_iter, tol, verbose)
  gamma <- .spca_scalar(gamma, "gamma")
  fit <- .spca_robust_solve(prepared, controls, gamma, penalty, B_init)
  out <- .spca_result(fit, prepared, "robspca", alpha * fit$norm2,
                      beta * fit$norm2, penalty)
  out$gamma <- gamma
  out
}
