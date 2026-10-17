#' @rdname spca_path
#' @export
spca_path <- function(X, k = NULL, alpha = c(1e-4, 1e-3, 1e-2),
                      beta = 1e-4, center = TRUE, scale = FALSE,
                      max_iter = 1000, tol = 1e-5,
                      method = c("spca", "rspca"), o = 20, q = 2,
                      verbose = TRUE, penalty = c("elasticnet", "l0"),
                      B_init = NULL) {
  method <- match.arg(method)
  penalty <- match.arg(penalty)
  if (!is.numeric(alpha) || is.complex(alpha) || !length(alpha) ||
      any(!is.finite(alpha)) || any(alpha < 0)) {
    stop("alpha must be a nonempty vector of finite nonnegative numbers.",
         call. = FALSE)
  }
  prepared <- .spca_preprocess(X, center, scale, sparse = method == "rspca")
  controls <- .spca_controls(k, prepared$n, prepared$p, alpha[1L], beta,
                              max_iter, tol, verbose)
  factor <- .spca_factor(prepared, controls$k, randomized = method == "rspca",
                         o = o, q = q)
  fits <- vector("list", length(alpha))
  nonzero <- matrix(0L, length(alpha), controls$k,
                    dimnames = list(NULL, paste0("PC", seq_len(controls$k))))
  selected_features <- iterations <- integer(length(alpha))
  for (i in seq_along(alpha)) {
    controls$alpha <- alpha[i]
    fit <- .spca_solve(factor, controls, penalty, B_init)
    fits[[i]] <- .spca_result(fit, prepared, method, alpha[i] * factor$norm2,
                              beta * factor$norm2, penalty)
    B_init <- fit$B
    nonzero[i, ] <- colSums(B_init != 0)
    selected_features[i] <- sum(rowSums(B_init != 0) > 0L)
    iterations[i] <- fits[[i]]$iterations
  }
  structure(list(fits = fits, alpha = alpha, nonzero = nonzero,
                  selected_features = selected_features, iterations = iterations,
                  method = method, penalty = penalty), class = "spca_path")
}

#' @rdname spca_path
#' @export
print.spca_path <- function(x, ...) {
  out <- data.frame(alpha = x$alpha, selected_features = x$selected_features,
                    iterations = x$iterations,
                    converged = vapply(x$fits, function(fit) fit$converged, logical(1)))
  print(out, row.names = FALSE)
  invisible(x)
}
