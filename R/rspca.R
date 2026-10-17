#' @rdname rspca
#' @export
rspca <- function(X, k = NULL, alpha = 1e-4, beta = 1e-4, center = TRUE, scale = FALSE, max_iter = 1000, tol = 1e-5, o = 20, q = 2, verbose = TRUE, penalty = c("elasticnet", "l0"), B_init = NULL) {
  UseMethod("rspca")
}

#' @export
rspca.default <- function(X, k = NULL, alpha = 1e-4, beta = 1e-4, center = TRUE, scale = FALSE, max_iter = 1000, tol = 1e-5, o = 20, q = 2, verbose = TRUE, penalty = c("elasticnet", "l0"), B_init = NULL) {
  penalty <- match.arg(penalty)
  prepared <- .spca_preprocess(X, center, scale, sparse = TRUE)
  controls <- .spca_controls(k, prepared$n, prepared$p, alpha, beta,
                              max_iter, tol, verbose)
  factor <- .spca_factor(prepared, controls$k, randomized = TRUE, o = o, q = q)
  fit <- .spca_solve(factor, controls, penalty, B_init)
  .spca_result(fit, prepared, "rspca", alpha * factor$norm2,
               beta * factor$norm2, penalty)
}
