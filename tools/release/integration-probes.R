# Exercise the sparsepca calls actually consumed by current reverse dependencies.
# This complements complete checks of the four CRAN strong reverse dependencies.
# It covers the optional parameters integration and both scPCA solver choices.
#
# Usage: Rscript tools/release/integration-probes.R OUTPUT_DIR SPARSEPCA_VERSION
# Select baseline/candidate with R_LIBS_USER before starting R; do not unload and
# swap namespaces within one process. All input data below are reproducible.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: integration-probes.R OUTPUT_DIR SPARSEPCA_VERSION")
}
output_dir <- args[1L]
expected_version <- args[2L]
stopifnot(as.character(utils::packageVersion("sparsepca")) == expected_version)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))

results <- list()
run_case <- function(name, expression) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  failed <- inherits(value, "error")
  results[[name]] <<- list(
    status = if (failed) "ERROR" else "OK", warnings = warnings,
    error = if (failed) conditionMessage(value) else NULL,
    value = if (failed) NULL else value)
  cat(name, results[[name]]$status, "\n")
  if (failed) cat(conditionMessage(value), "\n")
  if (length(warnings)) cat(paste(warnings, collapse = "\n"), "\n")
  invisible(value)
}

set.seed(761)
X <- matrix(rnorm(120L * 8L), 120L, 8L)
X[, 2L] <- 0.8 * X[, 1L] + 0.2 * X[, 2L]
X[, 4L] <- 0.7 * X[, 3L] + 0.3 * X[, 4L]
colnames(X) <- paste0("x", seq_len(ncol(X)))

run_case("categoryEncodings::encode_SPCA", {
  dat <- as.data.frame(X[, seq_len(5L), drop = FALSE])
  dat$category <- rep(letters[seq_len(12L)], each = 10L)
  ans <- categoryEncodings::encode_SPCA(dat, fact = "category", keep_factor = FALSE)
  stopifnot(nrow(ans) == nrow(dat), all(vapply(ans, is.numeric, logical(1))))
  ans
})

for (center in c(FALSE, TRUE)) {
  for (scale in c(FALSE, TRUE)) {
    run_case(paste0("clusterability scores center=", center, " scale=", scale), {
      ans <- getFromNamespace("perform_spca_sparsepca", "clusterability")(
        X, center = center, scale = scale, alpha = 0.01, beta = 0.001)
      stopifnot(length(ans) == nrow(X), all(is.finite(ans)))
      ans
    })
  }
}
run_case("clusterability::clusterabilitytest", {
  ans <- clusterability::clusterabilitytest(X, test = "dip", reduction = "spca",
                                          spca_method = "VP")
  stopifnot(is.list(ans))
  ans
})

run_case("robCompositions::spca_logrs", {
  positive <- exp(X[, seq_len(5L), drop = FALSE])
  ans <- robCompositions::spca_logrs(positive, alpha = 0.01, k = 2L, draw = FALSE)
  stopifnot(identical(dim(ans$model$loadings), c(20L, 2L)),
            identical(dim(ans$model$scores), c(120L, 2L)),
            is.matrix(ans[["model summary"]]),
            length(ans$expl.var) == 2L, all(is.finite(ans$expl.var)))
  ans
})

run_case("SparseBiplots::ElasticNet_HJBiplot", {
  ans <- SparseBiplots::ElasticNet_HJBiplot(
    as.data.frame(X[, seq_len(6L), drop = FALSE]), Lambda = 0.0001, Alpha = 0.001)
  stopifnot(identical(dim(ans$loadings), c(6L, 6L)))
  ans
})

run_case("parameters::principal_components sparse", {
  ans <- parameters::principal_components(as.data.frame(X), n = 2L, sparse = TRUE,
    k = 2L, alpha = 0.01, beta = 0.001, threshold = 0)
  model <- attr(ans, "model")
  stopifnot(is.data.frame(ans), length(model$sdev) == 2L,
            all(is.finite(model$sdev)), identical(dim(model$scores), c(120L, 2L)))
  ans
})

run_case("parameters::principal_components robust", {
  contaminated <- X
  contaminated[cbind(c(3L, 20L, 47L, 80L), c(1L, 3L, 6L, 8L))] <- 12
  ans <- parameters::principal_components(as.data.frame(contaminated),
    n = 2L, sparse = "robust", k = 2L, alpha = 0.01, beta = 0.001,
    gamma = 0.5, max_iter = 500L, threshold = 0)
  model <- attr(ans, "model")
  stopifnot(is.data.frame(ans), all(is.finite(model$sdev)),
            identical(dim(model$sparse), c(120L, 8L)))
  ans
})

for (algorithm in c("var_proj", "rand_var_proj")) {
  run_case(paste0("scPCA::scPCA ", algorithm), {
    set.seed(119)
    background <- matrix(rnorm(100L * 8L), 100L, 8L)
    colnames(background) <- colnames(X)
    ans <- scPCA::scPCA(target = X, background = background,
      n_eigen = 2L, alg = algorithm, contrasts = 0.1, penalties = 0.01,
      n_centers = 2L, parallel = FALSE)
    stopifnot(identical(dim(ans$rotation), c(8L, 2L)),
              identical(dim(ans$x), c(120L, 2L)), all(is.finite(ans$rotation)))
    ans
  })
}

saveRDS(results, file.path(output_dir, "probe-results.rds"))
summary <- do.call(rbind, lapply(names(results), function(name) {
  result <- results[[name]]
  data.frame(probe = name, sparsepca = expected_version, status = result$status,
             warnings = length(result$warnings),
             error = if (is.null(result$error)) "" else result$error,
             stringsAsFactors = FALSE)
}))
utils::write.csv(summary, file.path(output_dir, "probe-summary.csv"), row.names = FALSE)
packages <- c("sparsepca", "rsvd", "Matrix", "categoryEncodings", "clusterability",
              "robCompositions", "SparseBiplots", "parameters", "scPCA")
versions <- vapply(packages, function(package) {
  tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
}, character(1))
utils::write.csv(data.frame(package = packages, version = versions),
                 file.path(output_dir, "package-versions.csv"), row.names = FALSE)
if (any(summary$status == "ERROR")) quit(status = 1L)
