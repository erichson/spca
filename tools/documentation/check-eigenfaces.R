#!/usr/bin/env Rscript
# Recheck the saved application claims from the input and fitted factors.
# This does not refit PCA or SPCA. It reconstructs the randomized sketch once.
# Rscript tools/documentation/check-eigenfaces.R /path/to/faces.RData
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the path to faces.RData.")
if (!requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("rsvd", quietly = TRUE)) {
  stop("This verification requires digest and the recorded rsvd dependency.")
}
application <- readRDS("tools/documentation/results/eigenfaces-data.rds")
stricter_file <- "tools/documentation/results/eigenfaces-sparse-stricter.rds"
if (file.exists(stricter_file)) {
  stricter <- readRDS(stricter_file)
  stopifnot(identical(stricter$metadata$data_sha256, application$metadata$data_sha256),
              identical(stricter$metadata$runtime_md5, application$metadata$runtime_md5),
              identical(stricter$metadata$seed, application$metadata$seed),
              identical(stricter$parameters, application$settings$sparse),
              identical(stricter$objective[seq_along(application$fits$sparse$objective)],
                          application$fits$sparse$objective))
  # Add the separate completed sensitivity fit to verification only.
  # It remains distinct from the five fits used in the main application.
  application$fits$stricter_sparse <- stricter[c(
    "loadings", "transform", "scores", "center", "scale", "display_signs", "objective")]
  application$fit_metadata$stricter_sparse <- stricter$metadata
  application$settings$stricter_sparse <- stricter$parameters
  row <- stricter$metrics
  row$model <- "stricter_sparse"
  application$summary <- rbind(application$summary, row)
  B <- stricter$loadings
  application$component_metrics <- rbind(application$component_metrics,
    data.frame(model = "stricter_sparse", component = seq_len(ncol(B)),
                 nonzero_pixels = colSums(B != 0), nonzero_fraction = colMeans(B != 0),
                 loading_norm = sqrt(colSums(B^2)),
                 largest_absolute_loading = apply(abs(B), 2L, max),
                 display_sign = stricter$display_signs))
}
stopifnot(identical(digest::digest(file = args[1L], algo = "sha256"),
                      application$metadata$data_sha256))
env <- new.env(parent = emptyenv())
load(args[1L], envir = env)
X <- t(as.matrix(env$faces))
rm(env)
stopifnot(identical(dim(X), c(2410L, 8064L)))
reference <- application$fits$ordinary
processed <- scale(X, center = reference$center, scale = reference$scale)
rm(X)
energy <- sum(processed^2)
full_norm2 <- sum(reference$scores[, 1L]^2)
set.seed(application$metadata$seed)
sketch <- rsvd::rqb(processed, k = 45L, p = 0L, q = 2L)$B
sketch_norm2 <- svd(sketch, nu = 0L, nv = 0L)$d[1L]^2
checks <- list()
check <- function(model, metric, actual, expected, tolerance = 1e-9) {
  equal <- isTRUE(all.equal(actual, expected, tolerance = tolerance,
                             check.attributes = FALSE))
  checks[[length(checks) + 1L]] <<- data.frame(
    model = model, metric = metric, pass = equal,
    maximum_absolute_difference = max(abs(actual - expected)))
  if (!equal) stop("Verification failed: ", model, " / ", metric)
}
for (model in names(application$fits)) {
  fit <- application$fits[[model]]
  setting <- application$settings[[model]]
  target <- application$summary[application$summary$model == model, ]
  check(model, "orthogonal_transform", crossprod(fit$transform), diag(25L))
  check(model, "training_means", fit$center, reference$center)
  check(model, "training_scales", fit$scale, reference$scale)
  check(model, "nonzero_loadings", sum(fit$loadings != 0), target$nonzero_loadings, 0)
  check(model, "active_components", sum(colSums(fit$loadings^2) > 0), target$active_components, 0)
  coefficient_rows <- application$component_metrics[
    application$component_metrics$model == model, ]
  check(model, "component_support", colSums(fit$loadings != 0), coefficient_rows$nonzero_pixels, 0)
  check(model, "component_norms", sqrt(colSums(fit$loadings^2)), coefficient_rows$loading_norm)
  check(model, "largest_absolute_loading", apply(abs(fit$loadings), 2L, max),
          coefficient_rows$largest_absolute_loading)
  lowrank_sse <- complete_sse <- sparse_l1 <- sparse_sse <- sparse_count <- 0
  score_maximum_error <- 0
  for (start in seq.int(1L, nrow(processed), by = 100L)) {
    rows <- seq.int(start, min(nrow(processed), start + 99L))
    direct_scores <- processed[rows, , drop = FALSE] %*% fit$loadings
    score_maximum_error <- max(score_maximum_error,
                                max(abs(direct_scores - fit$scores[rows, , drop = FALSE])))
    residual <- processed[rows, , drop = FALSE] -
      tcrossprod(fit$scores[rows, , drop = FALSE], fit$transform)
    lowrank_sse <- lowrank_sse + sum(residual^2)
    if (model == "robust") {
      S <- sign(residual) * pmax(abs(residual) - setting$gamma, 0)
      sparse_count <- sparse_count + sum(S != 0)
      sparse_l1 <- sparse_l1 + sum(abs(S))
      sparse_sse <- sparse_sse + sum(S^2)
      complete_sse <- complete_sse + sum((residual - S)^2)
      if (start == 1L) {
        check(model, "selected_sparse_errors",
                S[application$metadata$selected_rows, , drop = FALSE],
                application$selected_sparse)
      }
    } else complete_sse <- complete_sse + sum(residual^2)
  }
  check(model, "scores_from_processed_input", score_maximum_error, 0)
  check(model, "full_data_relative_error", sqrt(lowrank_sse / energy),
          target$relative_reconstruction_error)
  check(model, "complete_model_relative_error", sqrt(complete_sse / energy),
          target$relative_complete_model_error)
  check(model, "sparse_count", sparse_count, target$nonzero_sparse_errors, 0)
  check(model, "sparse_l1", sparse_l1, target$sparse_l1)
  check(model, "sparse_sse", sparse_sse, target$sparse_sse)
  norm2 <- if (setting$method == "rspca") sketch_norm2 else full_norm2
  penalty <- norm2 * (setting$alpha * sum(abs(fit$loadings)) +
                        0.5 * setting$beta * sum(fit$loadings^2))
  if (model == "robust") penalty <- penalty + setting$gamma * sparse_l1
  check(model, "full_data_objective", 0.5 * complete_sse + penalty,
          target$full_data_objective)
  if (length(fit$objective)) {
    check(model, "iteration_count", length(fit$objective), target$iterations, 0)
    check(model, "stored_final_objective", tail(fit$objective, 1L),
            target$solver_final_objective)
    loss <- if (setting$method == "rspca") {
      0.5 * sum((sketch - tcrossprod(sketch %*% fit$loadings, fit$transform))^2)
    } else 0.5 * complete_sse
    check(model, "solver_objective_from_factors", loss + penalty,
            target$solver_final_objective)
    stopifnot(all(diff(fit$objective) <= 64 * .Machine$double.eps * energy))
    improvement <- (fit$objective[length(fit$objective) - 1L] - tail(fit$objective, 1L)) /
      tail(fit$objective, 1L)
    check(model, "relative_objective_improvement", improvement,
            target$final_relative_improvement)
    stopifnot(identical(improvement <= application$fit_metadata[[model]]$common$tol,
                         target$objective_converged))
  }
}
results <- do.call(rbind, checks)
print(results, row.names = FALSE)
cat("\nVerified", nrow(results), "claims from saved factors and the checksum-verified input.\n")
cat("R:", R.version.string, "\nrsvd:", as.character(utils::packageVersion("rsvd")), "\n")
