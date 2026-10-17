#!/usr/bin/env Rscript
# Generate the README's iris and sensor applications from public R APIs.
# From the repository root:
#   Rscript tools/documentation/application-examples.R
#   Rscript tools/documentation/application-examples.R --plot-only
# An optional final argument selects another repository root.
# Fitting uses the development dependency pkgload. Plot-only uses base R and
# the saved data/fit snapshot; it neither loads sparsepca nor refits models.
args <- commandArgs(trailingOnly = TRUE)
plot_only <- "--plot-only" %in% args
paths <- args[args != "--plot-only"]
if (length(paths) > 1L || any(grepl("^--", paths))) {
  stop("Usage: application-examples.R [--plot-only] [repository-root]")
}
root <- normalizePath(if (length(paths)) paths else ".")
result_dir <- file.path(root, "tools", "documentation", "results")
figure_dir <- file.path(root, "man", "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
snapshot <- file.path(result_dir, "application-data.rds")

# A sign convention is needed only to compare and label the pictures.
# Scores and coefficient columns must be flipped together. The fitted
# objects and the CSV columns named raw_* keep the original signs.
display_components <- function(loadings, scores) {
  signs <- vapply(seq_len(ncol(loadings)), function(j) {
    value <- loadings[which.max(abs(loadings[, j])), j]
    if (value < 0) -1 else 1
  }, numeric(1L))
  list(loadings = sweep(loadings, 2L, signs, "*"),
       scores = sweep(scores, 2L, signs, "*"), signs = signs)
}

if (plot_only) {
  if (!file.exists(snapshot)) stop("Run the fitting step before --plot-only.")
  application <- readRDS(snapshot)
} else {
  pkgload::load_all(root, quiet = TRUE)

  # Labels are deliberately kept outside both fitting calls.
  iris_X <- as.matrix(datasets::iris[, 1:4])
  iris_ordinary <- stats::prcomp(iris_X, rank. = 2,
                                 center = TRUE, scale. = TRUE)
  iris_sparse <- sparsepca::spca(iris_X, k = 2, alpha = 0.01, beta = 0.001,
                                 scale = TRUE, tol = 1e-8, max_iter = 3000,
                                 verbose = FALSE)
  iris_display <- list(
    ordinary = display_components(iris_ordinary$rotation, iris_ordinary$x),
    sparse = display_components(iris_sparse$loadings, iris_sparse$scores))
  iris_processed <- scale(iris_X, center = iris_sparse$center,
                            scale = iris_sparse$scale)
  iris_reconstruction <- list(
    ordinary = tcrossprod(iris_ordinary$x, iris_ordinary$rotation),
    sparse = tcrossprod(iris_sparse$scores, iris_sparse$transform))
  iris_error <- vapply(iris_reconstruction, function(estimate) {
    norm(estimate - iris_processed, "F") / norm(iris_processed, "F")
  }, numeric(1L))

  # Two shared signals are measured by separate groups of four sensors.
  # Independent small measurement errors are followed by four added spikes.
  set.seed(17)
  phase <- seq(0, 4 * pi, length.out = 120)
  latent <- cbind(sin(phase), cos(phase))
  weights <- matrix(0, 8, 2)
  weights[1:4, 1] <- 1
  weights[5:8, 2] <- 1
  truth <- tcrossprod(latent, weights)
  noise <- matrix(rnorm(120 * 8, sd = 0.08), 120, 8)
  observed <- truth + noise
  spots <- cbind(c(15, 45, 75, 105), c(1, 6, 3, 8))
  increments <- c(7, -7, 8, -8)
  added <- matrix(0, 120, 8)
  added[spots] <- increments
  observed[spots] <- observed[spots] + increments
  sensor_ordinary <- sparsepca::spca(observed, k = 2, alpha = 0.005,
                                     beta = 0.001, center = FALSE,
                                     tol = 1e-8, max_iter = 3000,
                                     verbose = FALSE)
  sensor_robust <- sparsepca::robspca(observed, k = 2, alpha = 0.005,
                                     beta = 0.001, gamma = 0.3,
                                     center = FALSE, tol = 1e-8,
                                     max_iter = 3000, verbose = FALSE)
  reconstruction <- list(
    ordinary = tcrossprod(sensor_ordinary$scores, sensor_ordinary$transform),
    robust = tcrossprod(sensor_robust$scores, sensor_robust$transform))
  relative_error <- vapply(reconstruction, function(estimate) {
    norm(estimate - truth, "F") / norm(truth, "F")
  }, numeric(1L))
  if (!iris_sparse$converged || !sensor_ordinary$converged ||
      !sensor_robust$converged) {
    stop("An application fit did not meet its objective stopping criterion.")
  }

  source_revision <- tryCatch(
    system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"),
              stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_)
  application <- list(
    format_version = 1L,
    metadata = list(
      generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      sparsepca_version = as.character(utils::packageVersion("sparsepca")),
      source_revision_at_generation = source_revision,
      R_version = R.version.string,
      RNG_kind = RNGkind(),
      thread_environment = Sys.getenv(c("OPENBLAS_NUM_THREADS",
                                          "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
      generator = "tools/documentation/application-examples.R"),
    iris = list(X = iris_X, species = datasets::iris$Species,
                ordinary = iris_ordinary, sparse = iris_sparse,
                display = iris_display,
                reconstruction = iris_reconstruction,
                relative_reconstruction_error = iris_error,
                parameters = list(k = 2, alpha = 0.01, beta = 0.001,
                                    center = TRUE, scale = TRUE,
                                    tol = 1e-8, max_iter = 3000)),
    sensors = list(phase = phase, latent = latent, weights = weights,
                   truth = truth, noise = noise, observed = observed,
                   added = added, spots = spots,
                   ordinary = sensor_ordinary, robust = sensor_robust,
                   reconstruction = reconstruction,
                   relative_error = relative_error,
                   parameters = list(seed = 17, n = 120, p = 8, k = 2,
                                       noise_sd = 0.08, alpha = 0.005,
                                       beta = 0.001, gamma = 0.3,
                                       center = FALSE, scale = FALSE,
                                       tol = 1e-8, max_iter = 3000)))
  saveRDS(application, snapshot, compress = "xz", version = 3)

  iris_scores <- do.call(rbind, lapply(c("ordinary", "sparse"), function(method) {
    raw <- if (method == "ordinary") iris_ordinary$x else iris_sparse$scores
    data.frame(method = method, flower = seq_len(nrow(iris_X)),
                 species = datasets::iris$Species,
                 raw_component_1 = raw[, 1], raw_component_2 = raw[, 2],
                 display_component_1 = iris_display[[method]]$scores[, 1],
                 display_component_2 = iris_display[[method]]$scores[, 2])
  }))
  iris_loadings <- do.call(rbind, lapply(c("ordinary", "sparse"), function(method) {
    raw <- if (method == "ordinary") iris_ordinary$rotation else iris_sparse$loadings
    data.frame(method = method, variable = rep(colnames(iris_X), 2),
                 component = rep(1:2, each = ncol(iris_X)),
                 raw_loading = as.vector(raw),
                 display_loading = as.vector(iris_display[[method]]$loadings))
  }))
  sensor_traces <- data.frame(
    sample = rep(seq_len(nrow(observed)), ncol(observed)),
    phase = rep(phase, ncol(observed)),
    sensor = rep(seq_len(ncol(observed)), each = nrow(observed)),
    truth = as.vector(truth), observed = as.vector(observed),
    ordinary_reconstruction = as.vector(reconstruction$ordinary),
    robust_reconstruction = as.vector(reconstruction$robust),
    added_error = as.vector(added),
    estimated_sparse_error = as.vector(sensor_robust$sparse))
  errors_to_show <- sensor_traces$added_error != 0 |
    sensor_traces$estimated_sparse_error != 0
  sensor_corruptions <- sensor_traces[errors_to_show,
    c("sample", "sensor", "added_error", "estimated_sparse_error")]

  metric_rows <- function(application, method, values) {
    data.frame(application = application, method = method,
                 metric = names(values), value = unname(values))
  }
  summary <- rbind(
    metric_rows("iris", "prcomp", c(
      nonzero_loadings = sum(iris_ordinary$rotation != 0),
      nonzero_component_1 = sum(iris_ordinary$rotation[, 1] != 0),
      nonzero_component_2 = sum(iris_ordinary$rotation[, 2] != 0),
      relative_reconstruction_error = unname(iris_error["ordinary"]))),
    metric_rows("iris", "spca", c(
      nonzero_loadings = sum(iris_sparse$loadings != 0),
      nonzero_component_1 = sum(iris_sparse$loadings[, 1] != 0),
      nonzero_component_2 = sum(iris_sparse$loadings[, 2] != 0),
      relative_reconstruction_error = unname(iris_error["sparse"]),
      iterations = iris_sparse$iterations,
      objective_converged = as.integer(iris_sparse$converged))),
    metric_rows("sensors", "spca", c(
      relative_clean_signal_error = relative_error["ordinary"],
      iterations = sensor_ordinary$iterations,
      objective_converged = as.integer(sensor_ordinary$converged))),
    metric_rows("sensors", "robspca", c(
      relative_clean_signal_error = relative_error["robust"],
      iterations = sensor_robust$iterations,
      objective_converged = as.integer(sensor_robust$converged),
      nonzero_sparse_errors = sum(sensor_robust$sparse != 0),
      injected_errors = sum(added != 0),
      injected_locations_with_nonzero_estimate =
        sum(sensor_robust$sparse[spots] != 0),
      other_nonzero_estimates = sum(sensor_robust$sparse[added == 0] != 0))))
  # Avoid inheriting the inner names of the two-element error vector.
  summary$metric <- sub("\\.(ordinary|robust)$", "", summary$metric)
  utils::write.csv(summary, file.path(result_dir, "application-summary.csv"),
                    row.names = FALSE)
  utils::write.csv(iris_scores, file.path(result_dir, "application-iris-scores.csv"),
                    row.names = FALSE)
  utils::write.csv(iris_loadings, file.path(result_dir, "application-iris-loadings.csv"),
                    row.names = FALSE)
  utils::write.csv(sensor_traces, file.path(result_dir, "application-sensor-traces.csv"),
                    row.names = FALSE)
  utils::write.csv(sensor_corruptions,
                    file.path(result_dir, "application-sensor-corruptions.csv"),
                    row.names = FALSE)
  session_lines <- c(capture.output(sessionInfo()), "", "Generation thread settings:",
                       paste(names(application$metadata$thread_environment),
                               application$metadata$thread_environment, sep = "="))
  writeLines(sub("[[:blank:]]+$", "", session_lines),
               file.path(result_dir, "application-session-info.txt"))
  print(summary, row.names = FALSE)
}

ink <- "#19384A"
muted <- "#526674"
grid_color <- "#E8EDF0"
ordinary_color <- "#BC6A32"
robust_color <- "#007B79"
truth_color <- "#233747"
species_colors <- c("#0072B2", "#D55E00", "#009E73")
species_shapes <- c(21, 22, 24)
coefficient_palette <- grDevices::colorRampPalette(
  c("#23638C", "#FFFFFF", "#B34949"))(201)

open_plot <- function(xlim, ylim, xlab = "", ylab = "") {
  plot(NA_real_, NA_real_, type = "n", xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, axes = FALSE, xaxs = "i", yaxs = "i")
  abline(h = pretty(ylim), v = pretty(xlim), col = grid_color, lwd = 0.7)
  axis(1, col = "#ADBBC4", col.axis = ink, col.ticks = "#ADBBC4")
  axis(2, las = 1, col = "#ADBBC4", col.axis = ink, col.ticks = "#ADBBC4")
  box(col = "#ADBBC4")
}

plot_iris <- function() {
  dat <- application$iris
  par(oma = c(0, 0, 4, 0), family = "sans", bg = "white",
        col.lab = ink, col.axis = ink, cex.axis = 0.86, cex.lab = 0.95)
  layout(matrix(c(1, 1, 2, 3, 4, 5, 6, 6), 4, 2, byrow = TRUE),
           heights = c(0.35, 3.5, 2.65, 0.68))
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = levels(dat$species), pch = species_shapes,
           pt.bg = species_colors, col = "white", pt.cex = 1.25,
           text.col = ink, bty = "n", horiz = TRUE, cex = 0.94)

  all_scores <- rbind(dat$display$ordinary$scores, dat$display$sparse$scores)
  score_xlim <- range(all_scores[, 1]) + c(-0.3, 0.3)
  score_ylim <- range(all_scores[, 2]) + c(-0.3, 0.3)
  for (method in c("ordinary", "sparse")) {
    par(mar = c(4.0, 4.0, 2.3, 1.15))
    scores <- dat$display[[method]]$scores
    open_plot(score_xlim, score_ylim, "Component 1 score", "Component 2 score")
    points(scores[, 1], scores[, 2], pch = species_shapes[dat$species],
             bg = grDevices::adjustcolor(species_colors[dat$species], alpha.f = 0.85),
             col = "white", lwd = 0.6, cex = 1.0)
    title(main = if (method == "ordinary") "Ordinary PCA" else "Sparse PCA",
            col.main = ink, cex.main = 1.1, line = 0.75)
  }

  for (method in c("ordinary", "sparse")) {
    par(mar = c(3.4, 6.7, 2.0, 1.15))
    B <- dat$display[[method]]$loadings
    image(1:2, 1:4, t(B[4:1, , drop = FALSE]), col = coefficient_palette,
            zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
            useRaster = FALSE)
    abline(h = (1:3) + 0.5, v = 1.5, col = "#DAE1E5", lwd = 0.7)
    for (j in 1:2) for (i in 1:4) {
      label <- if (B[i, j] == 0) "0" else sprintf("%.3f", B[i, j])
      text(j, 5 - i, label,
             col = if (abs(B[i, j]) > 0.6) "white" else ink,
             cex = 1.04, font = if (B[i, j] == 0) 2 else 1)
    }
    axis(1, at = 1:2, labels = c("Component 1", "Component 2"),
            tick = FALSE, cex.axis = 0.9)
    axis(2, at = 1:4, labels = rev(c("Sepal length", "Sepal width",
                                     "Petal length", "Petal width")),
            tick = FALSE, las = 1, cex.axis = 0.89)
    box(col = "#ADBBC4")
    title(main = "Coefficients used to compute the scores", line = 0.7,
            cex.main = 0.94, col.main = ink)
    mtext(sprintf("%d of 8 coefficients are nonzero", sum(B != 0)),
            side = 1, line = 2.25, cex = 0.89, col = muted)
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  edges <- seq(0.23, 0.77, length.out = length(coefficient_palette) + 1L)
  rect(edges[-length(edges)], 0.38, edges[-1], 0.64,
         col = coefficient_palette, border = NA)
  text(0.5, 0.9, "Coefficient applied to a standardized measurement",
         cex = 0.86, col = ink)
  text(c(0.23, 0.5, 0.77), 0.15, c("-1", "0", "+1"), cex = 0.83, col = muted)
  mtext("Iris: seeing groups, then reading the components", outer = TRUE,
          side = 3, line = 2.35, adj = 0.035, font = 2, cex = 1.32, col = ink)
  mtext("150 flowers; four standardized measurements; species labels are unused in fitting",
          outer = TRUE, side = 3, line = 0.7, adj = 0.035, cex = 0.89, col = muted)
}

plot_sensors <- function() {
  dat <- application$sensors
  index <- seq_len(nrow(dat$truth))
  par(oma = c(0.1, 0, 4.2, 0), family = "sans", bg = "white",
        col.lab = ink, col.axis = ink, cex.axis = 0.85, cex.lab = 0.91)
  layout(matrix(1:3, ncol = 1), heights = c(1.05, 1.1, 1.08))

  par(mar = c(2.3, 5.2, 2.45, 1.7))
  ylim <- range(dat$truth[, 1], dat$observed[, 1]) + c(-0.3, 1.1)
  open_plot(c(1, 120), ylim, ylab = "Sensor 1 reading")
  lines(index, dat$truth[, 1], col = truth_color, lwd = 1.6, lty = 2)
  points(index, dat$observed[, 1], pch = 16, cex = 0.55, col = "#80919C")
  points(15, dat$observed[15, 1], pch = 21, bg = ordinary_color,
           col = "white", cex = 1.35)
  text(18, dat$observed[15, 1], "Added spike: +7", adj = c(0, 0.5),
         cex = 0.86, col = ordinary_color)
  legend("topright", legend = c("Observed reading", "Known clean signal"),
           col = c("#80919C", truth_color), pch = c(16, NA),
           lty = c(NA, 2), lwd = c(NA, 1.6), bty = "n", cex = 0.84,
           text.col = ink)
  title(main = "An isolated bad reading interrupts a shared signal", line = 0.8,
          adj = 0, cex.main = 1.01, col.main = ink)

  par(mar = c(2.3, 5.2, 3.35, 1.7))
  ylim <- range(dat$truth[, 1], dat$reconstruction$ordinary[, 1],
                  dat$reconstruction$robust[, 1])
  ylim <- ylim + c(-0.15, 0.22) * diff(ylim)
  open_plot(c(1, 120), ylim, ylab = "Sensor 1 reconstruction")
  lines(index, dat$reconstruction$ordinary[, 1], col = ordinary_color, lwd = 1.5)
  lines(index, dat$reconstruction$robust[, 1], col = robust_color, lwd = 2.1)
  lines(index, dat$truth[, 1], col = truth_color, lwd = 1.2, lty = 2)
  legend("topright", legend = c("SPCA fit", "Robust SPCA fit", "Known clean signal"),
           col = c(ordinary_color, robust_color, truth_color), lty = c(1, 1, 2),
           lwd = c(1.5, 2.1, 1.2), bty = "n", cex = 0.82, text.col = ink)
  title(main = "The two fitted reconstructions (expanded vertical scale)",
          line = 1.7, adj = 0, cex.main = 1.01, col.main = ink)
  mtext(sprintf("Relative clean-signal error across all 8 sensors: SPCA %.1f%%; robust SPCA %.1f%%",
                  100 * dat$relative_error["ordinary"],
                  100 * dat$relative_error["robust"]),
          side = 3, line = 0.35, adj = 0, cex = 0.77, col = muted)

  par(mar = c(4.1, 5.2, 3.35, 1.7))
  error_limit <- max(abs(dat$added), abs(dat$robust$sparse))
  open_plot(c(1, 120), c(-1.2, 1.2) * error_limit,
              xlab = "Sample index", ylab = "Sparse error")
  abline(h = 0, col = "#A4B2BB", lwd = 1)
  estimated <- which(dat$robust$sparse != 0, arr.ind = TRUE)
  estimated_values <- dat$robust$sparse[estimated]
  segments(estimated[, 1], 0, estimated[, 1], estimated_values,
             col = robust_color, lwd = 1.5)
  points(estimated[, 1], estimated_values, pch = 16, col = robust_color, cex = 0.85)
  points(dat$spots[, 1], dat$added[dat$spots], pch = 1,
           col = truth_color, lwd = 1.3, cex = 1.15)
  text(dat$spots[, 1] + 2.5, 0.8 * dat$added[dat$spots],
         labels = paste("Sensor", dat$spots[, 2]),
         adj = c(0, 0.5), cex = 0.79, col = muted)
  legend("topright", legend = c("Added error", "Estimated S"),
           col = c(truth_color, robust_color), pch = c(1, 16),
           lwd = c(NA, 1.5), lty = c(NA, 1),
           bty = "n", cex = 0.82, text.col = ink)
  title(main = "The robust fit assigns large deviations to the sparse error matrix S",
          line = 1.7, adj = 0, cex.main = 1.01, col.main = ink)
  mtext(sprintf("All 8 sensors shown: %d entries of S are nonzero, including the 4 injected locations",
                  sum(dat$robust$sparse != 0)),
          side = 3, line = 0.35, adj = 0, cex = 0.77, col = muted)
  mtext("A few bad readings, a clearer shared signal", outer = TRUE,
          side = 3, line = 2.4, adj = 0.035, font = 2, cex = 1.31, col = ink)
  mtext("120 observations; 8 sensors; 2 latent signals; 4 injected spikes",
          outer = TRUE, side = 3, line = 0.75, adj = 0.035, cex = 0.9, col = muted)
}

write_figure <- function(name, plot_fun, width, height) {
  grDevices::svg(file.path(figure_dir, paste0(name, ".svg")),
                   width = width, height = height)
  tryCatch(plot_fun(), finally = grDevices::dev.off())
  grDevices::png(file.path(figure_dir, paste0(name, ".png")),
                   width = round(150 * width), height = round(150 * height), res = 150)
  tryCatch(plot_fun(), finally = grDevices::dev.off())
}

write_figure("iris-application", plot_iris, width = 10.6, height = 8.7)
write_figure("robust-sensors", plot_sensors, width = 10.6, height = 9.3)
