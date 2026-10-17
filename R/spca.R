#' @rdname spca
#' @export
spca <- function(X, k = NULL, alpha = 1e-4, beta = 1e-4, center = TRUE, scale = FALSE, max_iter = 1000, tol = 1e-5, verbose = TRUE, penalty = c("elasticnet", "l0"), B_init = NULL) {
  UseMethod("spca")
}

#' @export
spca.default <- function(X, k = NULL, alpha = 1e-4, beta = 1e-4, center = TRUE, scale = FALSE, max_iter = 1000, tol = 1e-5, verbose = TRUE, penalty = c("elasticnet", "l0"), B_init = NULL) {
  penalty <- match.arg(penalty)
  prepared <- .spca_preprocess(X, center, scale, sparse = FALSE)
  controls <- .spca_controls(k, prepared$n, prepared$p, alpha, beta,
                              max_iter, tol, verbose)
  factor <- .spca_factor(prepared, controls$k, randomized = FALSE)
  fit <- .spca_solve(factor, controls, penalty, B_init)
  .spca_result(fit, prepared, "spca", alpha * factor$norm2,
               beta * factor$norm2, penalty)
}
