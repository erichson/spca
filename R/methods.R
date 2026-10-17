.spca_print <- function(x) {
  cat("Standard deviations:\n")
  print(round(x$sdev, 3))
  cat("\nEigenvalues (adjusted score variances):\n")
  print(round(x$eigenvalues, 3))
  cat("\nSparse loadings:\n")
  loadings <- round(x$loadings, 3)
  print(loadings)
  invisible(loadings)
}

.spca_summary <- function(object) {
  variance <- object$sdev^2
  proportion <- if (object$var == 0) rep(0, length(variance)) else
    variance / object$var
  out <- rbind(variance, object$sdev, proportion, cumsum(proportion))
  rownames(out) <- c("Explained variance", "Standard deviations",
                     "Proportion of variance", "Cumulative proportion")
  colnames(out) <- paste0("PC", seq_along(variance))
  round(out, 3)
}

#' @rdname spca-methods
#' @export
print.spca <- function(x, ...) .spca_print(x)
#' @rdname spca-methods
#' @export
print.rspca <- function(x, ...) .spca_print(x)
#' @rdname spca-methods
#' @export
print.robspca <- function(x, ...) .spca_print(x)
#' @rdname spca-methods
#' @export
summary.spca <- function(object, ...) .spca_summary(object)
#' @rdname spca-methods
#' @export
summary.rspca <- function(object, ...) .spca_summary(object)
#' @rdname spca-methods
#' @export
summary.robspca <- function(object, ...) .spca_summary(object)

.spca_predict <- function(object, newdata) {
  X <- .spca_matrix(newdata, omit = FALSE, sparse = TRUE, min_rows = 0L)
  p <- nrow(object$loadings)
  if (ncol(X) != p) stop("newdata must have ", p, " columns.", call. = FALSE)
  training_names <- object$feature_names
  if (is.null(training_names) && is.numeric(object$center)) {
    training_names <- names(object$center)
  }
  if (!is.null(training_names) && !is.null(colnames(X))) {
    if (anyDuplicated(training_names) || anyDuplicated(colnames(X)) ||
        !setequal(training_names, colnames(X))) {
      stop("newdata must contain the same unique column names as the training data.",
           call. = FALSE)
    }
    X <- X[, match(training_names, colnames(X)), drop = FALSE]
  }
  if (inherits(X, "sparseMatrix")) {
    op <- .spca_sparse_products(X, object$center, object$scale)
    out <- op$forward(object$loadings)
  } else {
    if (!identical(object$center, FALSE)) X <- sweep(X, 2L, object$center, "-")
    if (!identical(object$scale, FALSE)) X <- sweep(X, 2L, object$scale, "/")
    out <- X %*% object$loadings
  }
  rownames(out) <- rownames(X)
  # Older Matrix releases can retain list(NULL, NULL) after as.matrix().
  # Use the same unnamed dense-matrix representation as base multiplication.
  if (is.null(rownames(out)) && is.null(colnames(out))) dimnames(out) <- NULL
  out
}

#' @rdname predict.spca
#' @export
predict.spca <- function(object, newdata, ...) {
  if (missing(newdata)) return(object$scores)
  .spca_predict(object, newdata)
}
#' @rdname predict.spca
#' @export
predict.rspca <- function(object, newdata, ...) {
  if (missing(newdata)) return(object$scores)
  .spca_predict(object, newdata)
}
#' @rdname predict.spca
#' @export
predict.robspca <- function(object, newdata, ...) {
  if (missing(newdata)) return(object$scores)
  .spca_predict(object, newdata)
}
