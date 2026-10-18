#' @rdname spca_target
#' @export
spca_target <- function(X, k = NULL, n_features, alpha = NULL,
                        beta = 1e-4, center = TRUE, scale = FALSE,
                        max_iter = 1000, tol = 1e-5,
                        method = c("spca", "rspca"), o = 20, q = 2,
                        verbose = TRUE, penalty = c("elasticnet", "l0"),
                        B_init = NULL, max_evals = 30L, count_tol = 0L) {
  method <- match.arg(method)
  penalty <- match.arg(penalty)
  n_features <- .spca_scalar(n_features, "n_features", 0, TRUE)
  count_tol <- .spca_scalar(count_tol, "count_tol", 0, TRUE)
  max_evals <- .spca_scalar(max_evals, "max_evals", 1, TRUE)
  if (max_evals > .Machine$integer.max || count_tol > .Machine$integer.max) {
    stop("max_evals and count_tol must not exceed the maximum integer.",
         call. = FALSE)
  }
  max_evals <- as.integer(max_evals)
  count_tol <- as.integer(count_tol)
  if (is.null(alpha)) {
    m <- min(max_evals, 13L)
    alpha <- if (m == 1L) 0 else if (m == 2L) c(0, 1) else
      c(0, 10^seq(-5, 0, length.out = m - 1L))
  } else if (!is.numeric(alpha) || is.complex(alpha) || !length(alpha) ||
             any(!is.finite(alpha)) || any(alpha < 0)) {
    stop("alpha must be a nonempty vector of finite nonnegative numbers.",
         call. = FALSE)
  }
  if (length(alpha) > max_evals) {
    stop("The initial alpha grid cannot contain more than max_evals values.",
         call. = FALSE)
  }
  alpha <- as.numeric(alpha)
  prepared <- .spca_preprocess(X, center, scale, sparse = method == "rspca")
  if (n_features > prepared$p) {
    stop("n_features cannot exceed the number of variables in X.", call. = FALSE)
  }
  n_features <- as.integer(n_features)
  controls <- .spca_controls(k, prepared$n, prepared$p, alpha[1L], beta,
                              max_iter, tol, verbose)
  factor <- .spca_factor(prepared, controls$k, randomized = method == "rspca",
                         o = o, q = q)

  tried <- numeric(max_evals)
  counts <- iterations <- warm_start <- integer(max_evals)
  converged <- logical(max_evals)
  nonzero <- matrix(0L, max_evals, controls$k,
                    dimnames = list(NULL, paste0("PC", seq_len(controls$k))))
  positive <- alpha[alpha > 0]
  reference <- if (length(positive)) min(positive) else 1
  best <- best_any <- NA_integer_
  best_fit <- best_any_fit <- NULL
  used <- 0L
  stop_reason <- "max_evals"

  for (i in seq_len(max_evals)) {
    if (i <= length(alpha)) {
      value <- alpha[i]
      warm_start[i] <- i - 1L
    } else {
      rows <- seq_len(used)
      value <- .spca_target_next(tried[rows], counts[rows], n_features, reference)
      if (is.na(value)) {
        stop_reason <- "no_refinement"
        break
      }
      # Refinements start from a retained fit. Diagnostics record that exact
      # evaluation, including when no fit has converged yet.
      if (!is.na(best)) {
        B_init <- best_fit$B
        warm_start[i] <- best
      } else {
        B_init <- best_any_fit$B
        warm_start[i] <- best_any
      }
    }
    controls$alpha <- value
    fit <- .spca_solve(factor, controls, penalty, B_init)
    B_init <- fit$B
    tried[i] <- value
    nonzero[i, ] <- colSums(fit$B != 0)
    counts[i] <- sum(rowSums(fit$B != 0) > 0L)
    iterations[i] <- length(fit$objective)
    converged[i] <- fit$converged
    used <- i
    rows <- seq_len(used)
    best <- .spca_target_pick(counts[rows], converged[rows], n_features)
    best_any <- .spca_target_pick(counts[rows], rep(TRUE, used), n_features)
    if (!is.na(best) && best == i) best_fit <- fit
    if (best_any == i) best_any_fit <- fit
    if (!is.na(best) && abs(counts[best] - n_features) <= count_tol) {
      stop_reason <- "target_met"
      break
    }
  }

  rows <- seq_len(used)
  search <- data.frame(alpha = tried[rows], selected_features = counts[rows],
                       difference = counts[rows] - n_features,
                       iterations = iterations[rows], converged = converged[rows],
                       warm_start = warm_start[rows])
  selected <- NA_integer_
  selected_alpha <- NA_real_
  result <- NULL
  if (is.na(best)) {
    stop_reason <- "no_converged_fit"
    warning("No candidate fit converged; increase max_iter or reconsider tol. ",
            "No fit was selected. Inspect the search diagnostics.", call. = FALSE)
  } else {
    selected <- counts[best]
    selected_alpha <- tried[best]
    # Scores and variance diagnostics are formed only for the selected fit.
    result <- .spca_result(best_fit, prepared, method,
                            selected_alpha * factor$norm2,
                            beta * factor$norm2, penalty)
  }
  structure(list(fit = result, n_features = n_features,
                  selected_features = selected, alpha = selected_alpha,
                  target_met = !is.na(best) &&
                    abs(selected - n_features) <= count_tol,
                  count_tol = count_tol, index = best, search = search,
                  nonzero = nonzero[rows, , drop = FALSE],
                  method = method, penalty = penalty, stop_reason = stop_reason),
            class = "spca_target")
}

# Compare support counts, never objectives from different penalties. On a
# distance tie, prefer fewer variables, then the earlier evaluation.
.spca_target_pick <- function(counts, converged, target) {
  eligible <- which(converged)
  if (!length(eligible)) return(NA_integer_)
  eligible[order(abs(counts[eligible] - target), counts[eligible], eligible)[1L]]
}

# Refine all observed intervals without assuming that support is monotone in
# alpha. Prioritize observed crossings of the target, then promising endpoint
# counts, splitting wider intervals first on ties. A fixed log1p scale gives
# the interval next to zero a finite width and reproducible ordering.
.spca_target_next <- function(alpha, counts, target, reference) {
  values <- sort(unique(alpha))
  if (length(values) < 2L) return(NA_real_)
  lower <- values[-length(values)]
  upper <- values[-1L]
  midpoint <- upper / 2
  positive <- lower > 0
  midpoint[positive] <- exp((log(lower[positive]) + log(upper[positive])) / 2)
  valid <- is.finite(midpoint) & midpoint > lower & midpoint < upper &
    !(midpoint %in% alpha)
  if (!any(valid)) return(NA_real_)
  distances <- vapply(values, function(value) {
    min(abs(counts[alpha == value] - target))
  }, numeric(1))
  smallest <- vapply(values, function(value) {
    min(counts[alpha == value])
  }, numeric(1))
  largest <- vapply(values, function(value) {
    max(counts[alpha == value])
  }, numeric(1))
  crossing <- pmin(smallest[-length(values)], smallest[-1L]) <= target &
    pmax(largest[-length(values)], largest[-1L]) >= target
  distance <- pmin(distances[-length(values)], distances[-1L])
  transformed <- numeric(length(values))
  large <- values > reference
  transformed[large] <- log(values[large]) - log(reference) +
    log1p(reference / values[large])
  transformed[!large] <- log1p(values[!large] / reference)
  width <- diff(transformed)
  candidates <- which(valid)
  chosen <- candidates[order(!crossing[candidates], distance[candidates],
                              -width[candidates],
                              lower[candidates])[1L]]
  midpoint[chosen]
}

#' @rdname spca_target
#' @export
print.spca_target <- function(x, ...) {
  cat("Sparse PCA feature targeting (", x$method, ", ", x$penalty, ")\n", sep = "")
  cat("Requested variables:", x$n_features,
      " | Count tolerance:", x$count_tol, "\n")
  if (is.null(x$fit)) {
    cat("No converged candidate; no fit selected.\n")
  } else {
    cat("Selected variables:", x$selected_features,
        " | Relative alpha:", format(x$alpha), "\n")
    cat(if (x$target_met) "Target met within count tolerance.\n" else
      "Closest converged fit found; target outside count tolerance.\n")
  }
  cat("Fits evaluated:", nrow(x$search), " | Stop reason:", x$stop_reason, "\n")
  invisible(x)
}
