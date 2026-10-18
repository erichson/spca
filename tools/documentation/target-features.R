#!/usr/bin/env Rscript
# Reproduce the README's three- and two-feature comparison.
# From the repository root:
#   Rscript tools/documentation/target-features.R
#   Rscript tools/documentation/target-features.R --plot-only
# An optional final argument selects another repository root. Fitting uses
# pkgload as a development tool. Plot-only loads saved results and uses base R.
args <- commandArgs(trailingOnly = TRUE)
plot_only <- "--plot-only" %in% args
paths <- args[args != "--plot-only"]
if (length(paths) > 1L || any(grepl("^--", paths))) {
  stop("Usage: target-features.R [--plot-only] [repository-root]")
}
root <- normalizePath(if (length(paths)) paths else ".")
result_dir <- file.path(root, "tools", "documentation", "results")
figure_dir <- file.path(root, "man", "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
snapshot <- file.path(result_dir, "target-features-iris.rds")

if (plot_only) {
  if (!file.exists(snapshot)) stop("Run the fitting step before --plot-only.")
  application <- readRDS(snapshot)
} else {
  pkgload::load_all(root, quiet = TRUE)
  X <- as.matrix(datasets::iris[, 1:4])
  fit_target <- function(requested) {
    parameters <- list(k = 2L, n_features = requested, beta = 0.001,
                     center = TRUE, scale = TRUE, tol = 1e-8,
                     max_iter = 3000L, method = "spca", max_evals = 30L,
                     count_tol = 0L, verbose = FALSE)
    target <- do.call(sparsepca::spca_target, c(list(X = X), parameters))
    if (is.null(target$fit) || !target$fit$converged || !target$target_met) {
      stop("The illustrative search did not produce a converged target fit.")
    }
    fit <- target$fit
    selected <- rowSums(fit$loadings != 0) > 0L
    if (sum(selected) != target$selected_features || sum(selected) != requested ||
        any(target$search$difference !=
              target$search$selected_features - target$n_features)) {
      stop("Stored feature counts disagree with the loading or search record.")
    }
    processed <- scale(X, center = fit$center, scale = fit$scale)
    reconstruction <- tcrossprod(fit$scores, fit$transform)
    relative_error <- norm(processed - reconstruction, "F") /
      norm(processed, "F")
    if (!isTRUE(all.equal(unname(fit$scores),
                         unname(processed %*% fit$loadings),
                         tolerance = 1e-12))) {
      stop("The selected scores disagree with the fitted loading map.")
    }
    # Preserve raw coefficients. Orient columns separately for legible display.
    signs <- vapply(seq_len(ncol(fit$loadings)), function(j) {
      value <- fit$loadings[which.max(abs(fit$loadings[, j])), j]
      if (value < 0) -1 else 1
    }, numeric(1L))
    display <- sweep(fit$loadings, 2L, signs, "*")
    list(parameters = parameters, target = target,
         selected_variable_names = colnames(X)[selected],
         display_loadings = display, display_signs = signs,
         relative_reconstruction_error = relative_error)
  }
  # Each request starts a separate search using the same deterministic settings.
  fits <- lapply(c(3L, 2L), fit_target)
  names(fits) <- c("three", "two")
  source_files <- c("R/spca_target.R", "R/solver.R", "R/utils.R",
                    "R/preprocess.R", "tools/documentation/target-features.R")
  source_md5 <- tools::md5sum(file.path(root, source_files))
  names(source_md5) <- source_files
  revision <- tryCatch(system2("git", c("-C", shQuote(root),
                                        "rev-parse", "HEAD"),
                                stdout = TRUE, stderr = FALSE),
                       error = function(e) NA_character_)
  git_state <- tryCatch(system2("git", c("-C", shQuote(root),
                                         "status", "--porcelain"),
                                 stdout = TRUE, stderr = FALSE),
                        error = function(e) "unavailable")
  application <- list(
    format_version = 2L,
    metadata = list(
      generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      sparsepca_version = as.character(utils::packageVersion("sparsepca")),
      source_revision_at_generation = revision,
      working_tree_at_generation = git_state,
      source_md5 = source_md5,
      R_version = R.version.string,
      thread_environment = Sys.getenv(c("OPENBLAS_NUM_THREADS",
                                         "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
      generator = "tools/documentation/target-features.R"),
    X = X, fits = fits)
  saveRDS(application, snapshot, compress = "xz", version = 3)
  summary <- do.call(rbind, lapply(fits, function(result) {
    target <- result$target
    fit <- target$fit
    data.frame(
      observations = nrow(X), variables = ncol(X), components = ncol(fit$loadings),
      requested_features = target$n_features,
      selected_features = target$selected_features,
      selected_variables = paste(result$selected_variable_names, collapse = ", "),
      target_met = target$target_met, count_tol = target$count_tol,
      selected_relative_alpha = target$alpha,
      evaluations = nrow(target$search), selected_evaluation = target$index,
      selected_fit_iterations = fit$iterations,
      selected_fit_converged = fit$converged,
      nonzero_coefficients = sum(fit$loadings != 0),
      relative_reconstruction_error = result$relative_reconstruction_error,
      stop_reason = target$stop_reason)
  }))
  searches <- do.call(rbind, lapply(fits, function(result) {
    data.frame(requested_features = result$target$n_features,
               result$target$search)
  }))
  coefficients <- do.call(rbind, lapply(fits, function(result) {
    B <- result$target$fit$loadings
    data.frame(
      requested_features = result$target$n_features,
      variable = rep(colnames(X), ncol(B)),
      component = rep(seq_len(ncol(B)), each = ncol(X)),
      raw_loading = as.vector(B),
      display_loading = as.vector(result$display_loadings))
  }))
  utils::write.csv(summary, file.path(result_dir, "target-features-summary.csv"),
                    row.names = FALSE)
  utils::write.csv(searches,
                    file.path(result_dir, "target-features-search.csv"),
                    row.names = FALSE)
  utils::write.csv(coefficients,
                    file.path(result_dir, "target-features-loadings.csv"),
                    row.names = FALSE)
  session <- c(capture.output(sessionInfo()), "", "Generation thread settings:",
                 paste(names(application$metadata$thread_environment),
                         application$metadata$thread_environment, sep = "="),
                 "", "Source revision at generation:", revision,
                 "Working tree at generation:", git_state,
                 "", "Source MD5 checksums (identify the working-tree code):",
                 paste(names(source_md5), unname(source_md5)))
  writeLines(sub("[[:blank:]]+$", "", session),
               file.path(result_dir, "target-features-session-info.txt"))
  print(summary, row.names = FALSE)
  print(searches, row.names = FALSE)
}

if (!identical(application$format_version, 2L)) {
  stop("Run the fitting step to save both feature-count searches.")
}

ink <- "#19384A"
muted <- "#526674"
grid_color <- "#E8EDF0"
fit_color <- "#23638C"
selected_color <- "#007B79"
palette <- grDevices::colorRampPalette(c("#23638C", "#FFFFFF", "#B34949"))(201)

draw_target <- function() {
  target <- application$fits$three$target
  search <- target$search
  B <- application$fits$three$display_loadings
  p <- nrow(B)
  n_eval <- nrow(search)
  par(oma = c(0, 0, 4.3, 0), family = "sans", bg = "white",
        col.lab = ink, col.axis = ink, cex.axis = 0.88, cex.lab = 0.95)
  layout(matrix(c(1, 2, 0, 3), 2, 2, byrow = TRUE),
           widths = c(1.05, 1.1), heights = c(1, 0.18))

  par(mar = c(4.4, 4.1, 2.8, 1.4))
  plot(NA_real_, NA_real_, type = "n", xlim = c(0.6, n_eval + 0.5),
         ylim = c(-0.15, p + 0.35), axes = FALSE,
         xlab = "Candidate evaluation", ylab = "Selected measurements",
         xaxs = "i", yaxs = "i")
  abline(h = 0:p, col = grid_color, lwd = 0.8)
  abline(h = target$n_features, col = selected_color, lty = 2, lwd = 1.4)
  axis(1, at = if (n_eval < 15L) seq_len(n_eval) else pretty(seq_len(n_eval)),
         col = "#ADBBC4", col.ticks = "#ADBBC4")
  axis(2, at = 0:p, las = 1, col = "#ADBBC4", col.ticks = "#ADBBC4")
  box(col = "#ADBBC4")
  lines(seq_len(n_eval), search$selected_features, col = fit_color, lwd = 1.3)
  points(seq_len(n_eval), search$selected_features, pch = 21,
           bg = ifelse(search$converged, fit_color, "white"),
           col = fit_color, cex = 1.0)
  points(target$index, target$selected_features, pch = 21,
           bg = selected_color, col = "white", cex = 1.55, lwd = 1.1)
  title(main = "Search for three measurements", line = 1.3,
          cex.main = 1.02, col.main = ink)
  mtext(sprintf("Selected relative alpha = %.5g", target$alpha),
          side = 3, line = 0.1, cex = 0.8, col = muted)
  legend("bottomleft", legend = c("Requested count", "Selected candidate",
                                    "Iteration limit reached"),
           col = c(selected_color, selected_color, fit_color),
           lty = c(2, NA, NA), pch = c(NA, 21, 21),
           pt.bg = c(selected_color, selected_color, "white"),
           lwd = c(1.4, NA, NA), pt.cex = 1.2,
           bty = "n", cex = 0.84, text.col = ink)

  par(mar = c(4.4, 7.1, 2.8, 1.0))
  limit <- max(1, abs(B))
  image(seq_len(ncol(B)), seq_len(p), t(B[p:1, , drop = FALSE]),
          col = palette, zlim = c(-limit, limit), axes = FALSE,
          xlab = "", ylab = "", useRaster = FALSE)
  abline(h = seq_len(p - 1L) + 0.5,
           v = seq_len(ncol(B) - 1L) + 0.5, col = "#DAE1E5", lwd = 0.7)
  for (j in seq_len(ncol(B))) for (i in seq_len(p)) {
    label <- if (B[i, j] == 0) "0" else if (abs(B[i, j]) < 0.005) {
      sprintf("%.4f", B[i, j])
    } else sprintf("%.3f", B[i, j])
    text(j, p + 1L - i, label,
           col = if (abs(B[i, j]) > 0.6 * limit) "white" else ink,
           cex = 1.03, font = if (B[i, j] == 0) 2 else 1)
  }
  axis(1, at = seq_len(ncol(B)), labels = paste("Component", seq_len(ncol(B))),
         tick = FALSE, cex.axis = 0.88)
  axis(2, at = seq_len(p), labels = rev(c("Sepal length", "Sepal width",
                                          "Petal length", "Petal width")),
         tick = FALSE, las = 1, cex.axis = 0.89)
  box(col = "#ADBBC4")
  title(main = "Three rows contribute to the scores", line = 1.3,
          cex.main = 1.02, col.main = ink)
  mtext(paste("Nonzero coefficients per component:",
                paste(colSums(B != 0), collapse = ", ")),
          side = 1, line = 2.2, cex = 0.82, col = muted)

  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  edges <- seq(0.15, 0.90, length.out = length(palette) + 1L)
  rect(edges[-length(edges)], 0.43, edges[-1], 0.67,
         col = palette, border = NA)
  text(0.525, 0.90, "Displayed loading coefficient", cex = 0.83, col = ink)
  text(c(0.15, 0.525, 0.90), 0.14,
         format(c(-limit, 0, limit), trim = TRUE), cex = 0.82, col = muted)
  mtext("Iris: choose how many measurements contribute", outer = TRUE,
          side = 3, line = 2.5, adj = 0.035, font = 2, cex = 1.27, col = ink)
  mtext("150 flowers; four standardized measurements; two score coordinates",
          outer = TRUE, side = 3, line = 0.85, adj = 0.035, cex = 0.9, col = muted)
}

grDevices::png(file.path(figure_dir, "target-features-iris.png"),
                 width = 1800, height = 900, res = 160, type = "cairo")
draw_target()
invisible(grDevices::dev.off())
grDevices::svg(file.path(figure_dir, "target-features-iris.svg"),
                 width = 11.25, height = 5.625, bg = "white", family = "sans")
draw_target()
invisible(grDevices::dev.off())
cat("Wrote target-features-iris.png and target-features-iris.svg.\n")
