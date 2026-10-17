args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare-inputs.R INPUT_ID OUTPUT_RDS")
input_id <- args[1]
output <- args[2]
set.seed(1729)

if (input_id %in% c("lowrank-small", "lowrank-medium", "lowrank-large")) {
  sizes <- list("lowrank-small" = c(2000L, 600L),
                "lowrank-medium" = c(5000L, 1500L),
                "lowrank-large" = c(10000L, 3000L))
  n <- sizes[[input_id]][1]
  p <- sizes[[input_id]][2]
  k <- 10L
  latent <- sweep(matrix(rnorm(n * k), n, k), 2, seq(3, 1, length.out = k), "*")
  loadings <- matrix(0, p, k)
  groups <- rep(seq_len(k), each = p / k)
  loadings[cbind(seq_len(p), groups)] <- 1
  X <- tcrossprod(latent, loadings) + matrix(rnorm(n * p, sd = 0.15), n, p)
} else if (input_id == "flat-small") {
  X <- matrix(rnorm(2000L * 600L), 2000L, 600L)
} else if (input_id == "sparse-large") {
  n <- 10000L
  p <- 3000L
  k <- 10L
  block <- p / k
  groups <- (seq_len(n) - 1L) %% k
  i <- rep(seq_len(n), each = block)
  j <- rep(groups * block, each = block) + rep(seq_len(block), times = n)
  x <- rep(runif(n, 0.8, 1.2), each = block) + rnorm(length(i), sd = 0.015)
  X <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, p))
} else {
  stop("Unknown input: ", input_id)
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(X, output, compress = FALSE)
cat(input_id, nrow(X), ncol(X), format(object.size(X), units = "MiB"), "\n")
