#!/usr/bin/env Rscript
# Reproduce the README illustration from self-contained synthetic data.
# Run from the repository root: Rscript tools/documentation/illustrate-loadings.R
# Development dependency: pkgload. No external data or images are required.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1L]) else normalizePath(".")
pkgload::load_all(root, quiet = TRUE)
result_dir <- file.path(root, "tools", "documentation", "results")
figure_dir <- file.path(root, "man", "figures")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20261017)
n <- 400L
p <- 12L
k <- 3L
latent <- matrix(rnorm(n * k), n, k)
latent[, 2L] <- 0.3 * latent[, 1L] + latent[, 2L]
latent[, 3L] <- 0.2 * latent[, 1L] + 0.3 * latent[, 2L] + latent[, 3L]
true_loadings <- kronecker(diag(k), matrix(1, 4L, 1L))
X <- tcrossprod(latent, true_loadings) + matrix(rnorm(n * p, sd = 0.12), n, p)
colnames(X) <- paste0(rep(c("A", "B", "C"), each = 4L), rep(1:4, k))
processed <- scale(X)
ordinary <- stats::prcomp(X, rank. = k, center = TRUE, scale. = TRUE)
elastic <- spca(X, k = k, alpha = 0.02, beta = 0.001,
                scale = TRUE, max_iter = 4000, tol = 1e-8, verbose = FALSE)
hard <- spca(X, k = k, alpha = 0.02, beta = 0.001, penalty = "l0",
             scale = TRUE, max_iter = 4000, tol = 1e-8, verbose = FALSE)
loadings <- list("Ordinary PCA" = ordinary$rotation,
                 "Elastic net SPCA" = elastic$loadings,
                 "Hard-threshold SPCA" = hard$loadings)
reconstructions <- list(tcrossprod(ordinary$x, ordinary$rotation),
                         tcrossprod(elastic$scores, elastic$transform),
                         tcrossprod(hard$scores, hard$transform))
errors <- vapply(reconstructions, function(reconstruction) {
  norm(processed - reconstruction, "F") / norm(processed, "F")
}, numeric(1L))
nonzero <- vapply(loadings, function(B) sum(B != 0), integer(1L))
summary <- data.frame(method = names(loadings), n = n, p = p, k = k,
                       nonzero = nonzero, total_coefficients = p * k,
                       relative_reconstruction_error = errors,
                       alpha = c(0, 0.02, 0.02), beta = c(0, 0.001, 0.001),
                       iterations = c(NA_integer_, elastic$iterations, hard$iterations),
                       converged = c(TRUE, elastic$converged, hard$converged))
utils::write.csv(summary, file.path(result_dir, "loadings-summary.csv"), row.names = FALSE)
coefficients <- do.call(rbind, lapply(seq_along(loadings), function(i) {
  data.frame(method = names(loadings)[i],
               variable = rep(colnames(X), k),
               component = rep(seq_len(k), each = p),
               loading = as.vector(loadings[[i]]))
}))
utils::write.csv(coefficients, file.path(result_dir, "loadings-coefficients.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "loadings-session-info.txt"))

plot_loadings <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  layout(matrix(c(1, 2, 3, 4, 4, 4), 2L, 3L, byrow = TRUE), heights = c(10, 2.2))
  par(oma = c(0, 0, 4.2, 0), family = "sans", bg = "white")
  bound <- max(abs(unlist(loadings)))
  palette <- grDevices::colorRampPalette(c("#23638C", "#FFFFFF", "#B34949"))(201)
  for (i in seq_along(loadings)) {
    par(mar = c(4.5, 3.2, 2.6, 1.1))
    B <- loadings[[i]]
    # Flip each component's sign for display only: largest coefficient positive.
    signs <- vapply(seq_len(ncol(B)), function(j) {
      sign(B[which.max(abs(B[, j])), j])
    }, numeric(1L))
    signs[signs == 0] <- 1
    B <- sweep(B, 2L, signs, "*")
    image(seq_len(k), seq_len(p), t(B[p:1L, , drop = FALSE]),
           zlim = c(-bound, bound), col = palette, axes = FALSE,
           xlab = "", ylab = "", useRaster = FALSE)
    axis(1, at = seq_len(k), labels = paste("Component", seq_len(k)),
          tick = FALSE, cex.axis = 0.83)
    axis(2, at = seq_len(p), labels = rev(colnames(X)),
          tick = FALSE, las = 1, cex.axis = 0.86)
    abline(h = c(4.5, 8.5), col = "#657482", lty = 3)
    box(col = "#BEC6CB")
    title(main = names(loadings)[i], cex.main = 1.0, col.main = "#173449")
    mtext(sprintf("%d / %d coefficients nonzero", nonzero[i], p * k),
            side = 1, line = 2.1, cex = 0.86, col = "#173449")
    mtext(sprintf("Reconstruction error: %.1f%%", 100 * errors[i]),
            side = 1, line = 3.25, cex = 0.8, col = "#536370")
  }
  par(mar = c(0.4, 0.6, 0.5, 0.6))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  edges <- seq(0.1, 0.9, length.out = length(palette) + 1L)
  rect(edges[-length(edges)], 0.35, edges[-1L], 0.67,
         col = palette, border = NA)
  text(0.5, 0.9, "Loading coefficient (common color scale)",
         cex = 0.85, col = "#173449")
  text(c(0.1, 0.5, 0.9), 0.1,
         c(sprintf("%.2f", -bound), "0", sprintf("%.2f", bound)), cex = 0.8)
  mtext("From shared variation to selected variables", outer = TRUE,
          side = 3, line = 2.3, adj = 0.03, cex = 1.32, font = 2, col = "#173449")
  mtext("400 observations; 12 variables in three correlated groups; k = 3",
          outer = TRUE, side = 3, line = 0.7, adj = 0.03, cex = 0.9, col = "#536370")
}

grDevices::svg(file.path(figure_dir, "sparse-loadings.svg"), width = 10.2, height = 7.0)
plot_loadings()
grDevices::dev.off()
grDevices::png(file.path(figure_dir, "sparse-loadings.png"), width = 1530, height = 1050, res = 150)
plot_loadings()
grDevices::dev.off()
print(summary, row.names = FALSE)
