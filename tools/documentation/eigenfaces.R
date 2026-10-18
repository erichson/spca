#!/usr/bin/env Rscript
# Reproduce the README face application through public fitting APIs.
# Data are downloaded separately; this script never fetches or redistributes them.
# See eigenfaces.md for the immutable data URL, checksum, and data provenance.
# Rscript tools/documentation/eigenfaces.R --fit faces.RData /path/to/checkpoints
# Rscript tools/documentation/eigenfaces.R --fit-model robust faces.RData /path/to/checkpoints
# Rscript tools/documentation/eigenfaces.R --assemble /path/to/checkpoints
# Rscript tools/documentation/eigenfaces.R --plot-only
# Fitting accepts --tol=VALUE for explicitly labeled sensitivity experiments.
# Run from the repository root. The plotting step uses saved derived data only.

args <- commandArgs(trailingOnly = TRUE)
tol_option <- grep("^--tol=", args)
if (length(tol_option) > 1L) stop("Supply at most one --tol=VALUE option.")
fit_tolerance <- if (length(tol_option))
  as.numeric(sub("^--tol=", "", args[tol_option])) else 1e-4
if (!is.finite(fit_tolerance) || fit_tolerance < 0) stop("Invalid fitting tolerance.")
if (length(tol_option)) args <- args[-tol_option]
mode <- if (length(args)) args[1L] else ""
if (!mode %in% c("--fit", "--fit-model", "--assemble", "--plot-only")) {
  stop("See the usage examples at the top of tools/documentation/eigenfaces.R")
}
if (length(tol_option) && !mode %in% c("--fit", "--fit-model")) {
  stop("--tol=VALUE applies only to fitting, not to saved results.")
}
root <- normalizePath(".")
result_dir <- file.path(root, "tools", "documentation", "results")
figure_dir <- file.path(root, "man", "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
snapshot <- file.path(result_dir, "eigenfaces-data.rds")
models <- c("ordinary", "sparse", "stronger", "robust", "robust_control")
settings <- list(
  ordinary = list(method = "prcomp", alpha = 0, beta = 0),
  sparse = list(method = "rspca", alpha = 1e-4, beta = 0.1),
  stronger = list(method = "rspca", alpha = 2e-4, beta = 0.2),
  robust = list(method = "robspca", alpha = 1e-4, beta = 0.01, gamma = 0.9),
  robust_control = list(method = "spca", alpha = 1e-4, beta = 0.01))
common <- list(k = 25L, center = TRUE, scale = TRUE,
                 tol = fit_tolerance, max_iter = 1000L)
source_sha256 <- "82674d0785f38c4c57eb688992eaf0fb2453e8c0409ef1f25b40da8de6a979a5"
data_url <- paste0("https://raw.githubusercontent.com/erichson/data/",
                    "fca24620977d21f5230ca9cc7a6107b25847d6fa/R/faces.RData")

if (mode %in% c("--fit", "--fit-model")) {
  expected_length <- if (mode == "--fit") 3L else 4L
  if (length(args) != expected_length) stop("Wrong number of fitting arguments.")
  selected <- if (mode == "--fit") models else args[2L]
  if (any(!selected %in% models)) stop("Unknown model name.")
  data_file <- normalizePath(args[expected_length - 1L])
  checkpoint_dir <- args[expected_length]
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- normalizePath(checkpoint_dir)
  if (!requireNamespace("digest", quietly = TRUE) ||
      !requireNamespace("pkgload", quietly = TRUE)) {
    stop("The fitting script requires the development packages digest and pkgload.")
  }
  actual_hash <- digest::digest(file = data_file, algo = "sha256")
  if (!identical(actual_hash, source_sha256)) stop("The face data checksum differs.")
  data_env <- new.env(parent = emptyenv())
  loaded <- load(data_file, envir = data_env)
  if (!identical(loaded, "faces") ||
      !identical(dim(data_env$faces), c(8064L, 2410L))) {
    stop("Unexpected face data objects or dimensions.")
  }
  X <- t(as.matrix(data_env$faces))
  rm(data_env)
  stopifnot(all(is.finite(X)), range(X) == c(0, 255))
  pkgload::load_all(root, quiet = TRUE)
  processed <- scale(X, center = TRUE, scale = TRUE)
  input_energy <- sum(processed^2)
  seed <- 42L
  metadata <- list(
    data_url = data_url, data_sha256 = source_sha256,
    data_git_blob = "884cf6fa86db4f26adbeb22eeb019101f99585c8",
    n = nrow(X), p = ncol(X), height = 96L, width = 84L,
    selected_rows = seq_len(25L), seed = seed, common = common,
    input_energy = input_energy,
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source_revision_at_generation = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    runtime_md5 = tools::md5sum(file.path(root, "R", c(
      "methods.R", "preprocess.R", "robspca.R", "rspca.R", "solver.R",
      "spca.R", "spca_path.R", "utils.R"))),
    sparsepca_version = as.character(utils::packageVersion("sparsepca")),
    R_version = R.version.string, RNG_kind = RNGkind(),
    thread_environment = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
    session = sub("[[:blank:]]+$", "", capture.output(sessionInfo())))

  for (model in selected) {
    params <- settings[[model]]
    message("Fitting ", model, " on ", nrow(X), " x ", ncol(X), " pixels")
    set.seed(seed)
    before <- proc.time()[["elapsed"]]
    if (model == "ordinary") {
      fit <- stats::prcomp(X, rank. = common$k, center = TRUE, scale. = TRUE)
      B <- A <- fit$rotation
      scores <- fit$x
      objective <- numeric()
      iterations <- NA_integer_
      converged <- NA
      improvement <- NA_real_
      effective_alpha <- effective_beta <- 0
    } else {
      call <- c(list(X = X), common, params[names(params) != "method"],
                  list(verbose = TRUE))
      if (params$method == "rspca") call <- c(call, list(o = 20L, q = 2L))
      fit <- do.call(getExportedValue("sparsepca", params$method), call)
      B <- fit$loadings
      A <- fit$transform
      scores <- fit$scores
      objective <- fit$objective
      iterations <- fit$iterations
      converged <- fit$converged
      improvement <- fit$improvement
      effective_alpha <- fit$alpha
      effective_beta <- fit$beta
    }
    elapsed <- proc.time()[["elapsed"]] - before
    # Measure the same full-data error for all methods, in bounded row blocks.
    # The robust model is also evaluated after including its fitted sparse term.
    reconstruction_sse <- complete_sse <- sparse_l1 <- sparse_sse <- 0
    nonzero_errors <- 0L
    for (start in seq.int(1L, nrow(X), by = 100L)) {
      rows <- seq.int(start, min(nrow(X), start + 99L))
      residual <- processed[rows, , drop = FALSE] -
        tcrossprod(scores[rows, , drop = FALSE], A)
      reconstruction_sse <- reconstruction_sse + sum(residual^2)
      if (model == "robust") {
        S <- fit$sparse[rows, , drop = FALSE]
        complete_sse <- complete_sse + sum((residual - S)^2)
        sparse_l1 <- sparse_l1 + sum(abs(S))
        sparse_sse <- sparse_sse + sum(S^2)
        nonzero_errors <- nonzero_errors + sum(S != 0)
        stopifnot(isTRUE(all.equal(S, sign(residual) * pmax(abs(residual) - params$gamma, 0),
                                    tolerance = 1e-12, check.attributes = FALSE)))
      } else complete_sse <- complete_sse + sum(residual^2)
    }
    full_objective <- 0.5 * complete_sse +
      effective_alpha * sum(abs(B)) + 0.5 * effective_beta * sum(B^2) +
      if (model == "robust") params$gamma * sparse_l1 else 0
    if (model %in% c("robust", "robust_control")) {
      stopifnot(isTRUE(all.equal(tail(objective, 1L), full_objective,
                                  tolerance = 1e-10, check.attributes = FALSE)))
    }
    if (length(objective) > 1L) {
      stopifnot(all(diff(objective) <= 64 * .Machine$double.eps * input_energy))
    }
    signs <- vapply(seq_len(ncol(B)), function(j) {
      if (B[which.max(abs(B[, j])), j] < 0) -1 else 1
    }, numeric(1L))
    metrics <- data.frame(
      model = model, method = params$method, n = nrow(X), p = ncol(X), k = common$k,
      alpha = params$alpha, beta = params$beta,
      gamma = if (model == "robust") params$gamma else NA_real_,
      nonzero_loadings = sum(B != 0), total_loadings = length(B),
      nonzero_loading_fraction = sum(B != 0) / length(B),
      active_components = sum(colSums(B^2) > 0),
      relative_reconstruction_error = sqrt(reconstruction_sse / input_energy),
      relative_complete_model_error = sqrt(complete_sse / input_energy),
      nonzero_sparse_errors = nonzero_errors,
      nonzero_sparse_error_fraction = nonzero_errors / length(processed),
      sparse_l1 = sparse_l1, sparse_sse = sparse_sse,
      full_data_objective = full_objective,
      solver_final_objective = if (length(objective)) tail(objective, 1L) else NA_real_,
      iterations = iterations, objective_converged = converged,
      final_relative_improvement = improvement, elapsed_seconds = elapsed)
    record <- list(metadata = metadata, parameters = params, metrics = metrics,
                     loadings = B, transform = A, scores = scores,
                     center = attr(processed, "scaled:center"),
                     scale = attr(processed, "scaled:scale"),
                     display_signs = signs, objective = objective,
                     selected_sparse = if (model == "robust")
                       fit$sparse[metadata$selected_rows, , drop = FALSE] else NULL)
    saveRDS(record, file.path(checkpoint_dir, paste0(model, ".rds")),
              compress = "gzip", version = 3)
    utils::write.csv(metrics, file.path(checkpoint_dir, paste0(model, "-summary.csv")),
                       row.names = FALSE)
    print(metrics, row.names = FALSE)
    rm(fit, B, A, scores, record)
    invisible(gc())
  }
  if (mode == "--fit-model") quit(status = 0L)
}

if (mode == "--assemble") {
  if (length(args) != 2L) stop("--assemble requires a checkpoint directory.")
  checkpoint_dir <- normalizePath(args[2L])
}
if (mode %in% c("--fit", "--assemble")) {
  records <- setNames(lapply(models, function(model) {
    readRDS(file.path(checkpoint_dir, paste0(model, ".rds")))
  }), models)
  for (record in records) {
    controls <- if (record$parameters$method == "prcomp")
      c("k", "center", "scale") else names(common)
    stopifnot(identical(record$metadata$data_sha256, source_sha256),
              identical(record$metadata$common[controls], common[controls]),
              identical(record$metadata$runtime_md5, records[[1L]]$metadata$runtime_md5))
  }
  summary <- do.call(rbind, lapply(records, `[[`, "metrics"))
  component_metrics <- do.call(rbind, lapply(models, function(model) {
    B <- records[[model]]$loadings
    data.frame(model = model, component = seq_len(ncol(B)),
                 nonzero_pixels = colSums(B != 0),
                 nonzero_fraction = colMeans(B != 0),
                 loading_norm = sqrt(colSums(B^2)),
                 largest_absolute_loading = apply(abs(B), 2L, max),
                 display_sign = records[[model]]$display_signs)
  }))
  objective_rows <- do.call(rbind, lapply(models, function(model) {
    obj <- records[[model]]$objective
    if (!length(obj)) return(NULL)
    data.frame(model = model, iteration = seq_along(obj), objective = obj)
  }))
  # Keep the fitted factors and selected errors so both plots and numerical
  # summaries can be inspected without repeating the expensive optimizations.
  # No source photographs or full dense sparse-error matrix are redistributed.
  application <- list(format_version = 2L, metadata = records$sparse$metadata,
                        fit_metadata = lapply(records, `[[`, "metadata"),
                        settings = settings, summary = summary,
                        component_metrics = component_metrics,
                        objective = objective_rows,
                        fits = lapply(records, function(record) {
                          record[c("loadings", "transform", "scores", "center",
                                   "scale", "display_signs", "objective")]
                        }),
                        selected_sparse = records$robust$selected_sparse)
  saveRDS(application, snapshot, compress = "xz", version = 3)
  utils::write.csv(summary, file.path(result_dir, "eigenfaces-summary.csv"), row.names = FALSE)
  utils::write.csv(component_metrics, file.path(result_dir, "eigenfaces-components.csv"), row.names = FALSE)
  utils::write.csv(objective_rows, file.path(result_dir, "eigenfaces-objectives.csv"), row.names = FALSE)
  session_lines <- c(application$metadata$session, "", "Thread settings:",
                     paste(names(application$metadata$thread_environment),
                             application$metadata$thread_environment, sep = "="))
  writeLines(session_lines, file.path(result_dir, "eigenfaces-session-info.txt"))
} else {
  if (length(args) != 1L) stop("--plot-only takes no additional arguments.")
  application <- readRDS(snapshot)
}

ink <- "#19384A"
muted <- "#526674"
palette <- grDevices::colorRampPalette(c("#235C85", "#FFFFFF", "#A93636"))(257)
display_loadings <- lapply(application$fits[1:4], function(fit) {
  magnitudes <- apply(abs(fit$loadings), 2L, max)
  sweep(sweep(fit$loadings, 2L, fit$display_signs, "*"),
           2L, ifelse(magnitudes > 0, magnitudes, 1), "/")
})
loading_limit <- 1
error_limit <- max(abs(application$selected_sparse))
display_order <- seq_len(25L)

plot_grid <- function(values, title, subtitle, labels, limit, legend_label) {
  layout(rbind(matrix(seq_len(25L), 5L, 5L, byrow = TRUE), rep(26L, 5L)),
           heights = c(rep(1, 5), 0.48))
  par(oma = c(0.5, 0.8, 4.4, 0.8), family = "sans", fg = ink, col = ink)
  for (j in display_order) {
    par(mar = c(0.20, 0.15, 1.45, 0.15))
    # The source vectors store 84 horizontal pixels before moving down a row.
    # Reversing the second index makes the top of the face the top of the plot.
    img <- matrix(values[, j], nrow = 84L, ncol = 96L)
    image(seq_len(84L), seq_len(96L), img[, 96L:1L],
            col = palette, zlim = c(-limit, limit), useRaster = TRUE,
            asp = 1, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    box(col = "#DCE4E9", lwd = 0.5)
    mtext(labels[j], side = 3, line = 0.3, cex = 0.72, col = muted)
  }
  par(mar = c(3.0, 8.5, 0.45, 8.5))
  image(seq(-limit, limit, length.out = length(palette)), c(0, 1),
          matrix(rep(seq(-limit, limit, length.out = length(palette)), 2), ncol = 2),
          col = palette, zlim = c(-limit, limit), useRaster = TRUE,
          axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  axis(1, at = c(-limit, 0, limit), labels = format(c(-limit, 0, limit), digits = 3),
         col = "#ADBBC4", col.axis = muted, cex.axis = 0.77, tck = -0.18)
  mtext(legend_label, side = 1, line = 2.1, cex = 0.73, col = muted)
  mtext(title, side = 3, outer = TRUE, line = 2.8, adj = 0.01,
          font = 2, cex = 1.3, col = ink)
  mtext(subtitle, side = 3, outer = TRUE, line = 1.2, adj = 0.01,
          cex = 0.87, col = muted)
}

write_figure <- function(name, plot_fun) {
  grDevices::svg(file.path(figure_dir, paste0(name, ".svg")), width = 10, height = 11.9)
  tryCatch(plot_fun(), finally = grDevices::dev.off())
  grDevices::png(file.path(figure_dir, paste0(name, ".png")),
                   width = 1600L, height = 1904L, res = 160)
  tryCatch(plot_fun(), finally = grDevices::dev.off())
}

titles <- c(ordinary = "Ordinary components use pixels throughout the face",
              sparse = "Sparse coefficients concentrate on smaller regions",
              stronger = "Stronger penalties reduce the number of active pixels",
              robust = "Robust fitting assigns large residuals to a separate term")
for (model in names(display_loadings)) {
  B <- application$fits[[model]]$loadings
  displayed <- display_loadings[[model]]
  params <- application$settings[[model]]
  metric <- application$summary[application$summary$model == model, ]
  parameter_text <- if (model == "ordinary") "PCA; 25 components" else
    paste0("alpha = ", format(params$alpha, scientific = TRUE),
             "; beta = ", params$beta,
             if (model == "robust") paste0("; gamma = ", params$gamma) else "")
  subtitle <- sprintf("%s; %.1f%% of coefficients nonzero", parameter_text,
                        100 * metric$nonzero_loading_fraction)
  labels <- sprintf("%02d   %.1f%% active", seq_len(ncol(B)), 100 * colMeans(B != 0))
  write_figure(paste0("eigenfaces-", model), function() {
    plot_grid(displayed, titles[model], subtitle, labels, loading_limit,
                "Coefficient relative to the largest magnitude in that component")
  })
}
S <- application$selected_sparse
labels <- sprintf("Image %02d   %.1f%% active", application$metadata$selected_rows,
                    100 * rowMeans(S != 0))
write_figure("eigenfaces-errors", function() {
  plot_grid(t(S), "The sparse term contains localized residual patterns",
              "First 25 images; fitted additive errors in standardized pixel units",
              labels, error_limit, "Estimated sparse error S; one common color scale across images")
})
