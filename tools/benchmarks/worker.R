# Run one fit in a fresh process. Input generation, loading, package loading,
# and warm-up are outside elapsed time. The RSS checkpoint is taken immediately
# after the fit, before blocked quality validation can change the high-water mark.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("Usage: worker.R SOURCE_DIR INPUT_RDS METHOD CASE REPETITION OUTPUT_CSV FINGERPRINT ENV_TXT")
}
source_dir <- args[1]
input_path <- args[2]
method <- args[3]
case <- args[4]
repetition <- as.integer(args[5])
output <- args[6]
fingerprint <- args[7]
environment_file <- args[8]
for (file in list.files(source_dir, pattern = "[.]R$", full.names = TRUE)) {
  sys.source(file, envir = .GlobalEnv)
}
invisible(loadNamespace("rsvd"))
invisible(loadNamespace("Matrix"))

pca <- startsWith(case, "pca-")
alpha <- if (pca) 0 else 0.001
beta <- if (pca) 0 else 0.0001
k <- 10L
o <- 20L
q <- 2L
max_iter <- 500L
tol <- 1e-5
fit_one <- function(X, warm = FALSE) {
  rank <- if (warm) 3L else k
  if (method == "prcomp") {
    return(stats::prcomp(X, rank. = rank, center = TRUE, scale. = FALSE))
  }
  if (method == "rspca-dense") X <- as.matrix(X)
  fun <- if (method == "spca") spca else rspca
  controls <- list(X = X, k = rank, alpha = alpha, beta = beta,
                   center = TRUE, scale = FALSE,
                   max_iter = if (warm) 5L else max_iter,
                   tol = tol, verbose = FALSE)
  if (method != "spca") controls <- c(controls, list(o = o, q = q))
  do.call(fun, controls)
}
set.seed(912)
warmup <- matrix(rnorm(80 * 30), 80, 30)
invisible(fit_one(warmup, warm = TRUE))
rm(warmup)
X <- readRDS(input_path)
invisible(gc())
set.seed(5000L + repetition)
warnings_seen <- character()
start <- proc.time()[["elapsed"]]
fit <- withCallingHandlers(fit_one(X), warning = function(w) {
  warnings_seen <<- c(warnings_seen, conditionMessage(w))
})
elapsed <- proc.time()[["elapsed"]] - start
status <- readLines("/proc/self/status")
rss_line <- status[startsWith(status, "VmHWM:")]
peak_rss_kib <- as.numeric(gsub("[^0-9]", "", rss_line))

if (method == "prcomp") {
  B <- fit$rotation
  A <- fit$rotation
  scores <- fit$x
  offset <- fit$center
  eigenvalues <- fit$sdev[seq_len(k)]^2
  norm2 <- fit$sdev[1]^2 * (nrow(X) - 1)
  alpha_effective <- beta_effective <- 0
  iterations <- 0L
  converged <- TRUE
  logged_objective <- NA_real_
} else {
  B <- fit$loadings
  A <- fit$transform
  scores <- fit$scores
  offset <- fit$center
  eigenvalues <- fit$eigenvalues
  norm2 <- if (alpha != 0) fit$alpha / alpha else NA_real_
  alpha_effective <- fit$alpha
  beta_effective <- fit$beta
  iterations <- fit$iterations
  converged <- fit$converged
  logged_objective <- tail(fit$objective, 1)
}

# Never densify the entire sparse input just to measure quality.
residual2 <- energy <- 0
block <- max(1L, floor(1e6 / ncol(X)))
for (first in seq.int(1L, nrow(X), by = block)) {
  rows <- seq.int(first, min(nrow(X), first + block - 1L))
  current <- as.matrix(X[rows, , drop = FALSE])
  if (!identical(offset, FALSE)) current <- sweep(current, 2L, offset, "-")
  residual <- current - tcrossprod(scores[rows, , drop = FALSE], A)
  residual2 <- residual2 + sum(residual^2)
  energy <- energy + sum(current^2)
}
nnz <- sum(B != 0)
l1 <- sum(abs(B))
l2_squared <- sum(B^2)
native_full_objective <- 0.5 * residual2 + alpha_effective * l1 +
  0.5 * beta_effective * l2_squared
record <- data.frame(
  case = case, method = method, repetition = repetition,
  n = nrow(X), p = ncol(X), k = k,
  input_class = class(X)[1], input_bytes = as.numeric(object.size(X)),
  input_density = if (inherits(X, "sparseMatrix")) length(X@x)/length(X) else 1,
  alpha = alpha, beta = beta, o = if (method == "spca" || method == "prcomp") NA else o,
  q = if (method == "spca" || method == "prcomp") NA else q,
  max_iter = max_iter, tol = tol,
  input_seed = 1729L, fit_seed = 5000L + repetition,
  elapsed_seconds = elapsed, peak_rss_mib = peak_rss_kib / 1024,
  iterations = iterations, converged = converged,
  reconstruction_relative_frobenius = sqrt(residual2 / energy),
  reconstruction_error_squared = residual2, input_energy = energy,
  fitted_energy_fraction = sum(scores^2) / energy,
  adjusted_score_energy_fraction = sum(eigenvalues) * (nrow(X) - 1) / energy,
  nonzeros = nnz, loading_zero_fraction = 1 - nnz/length(B),
  loading_l1 = l1, loading_l2_squared = l2_squared,
  spectral_norm_squared = norm2,
  alpha_effective = alpha_effective, beta_effective = beta_effective,
  logged_objective = logged_objective, full_objective_native = native_full_objective,
  warnings = length(warnings_seen), warning_text = paste(warnings_seen, collapse = " | "),
  source_fingerprint = fingerprint,
  stringsAsFactors = FALSE
)
write.csv(record, output, row.names = FALSE)
if (!file.exists(environment_file)) {
  capture.output({
    print(sessionInfo())
    cat("\nNumerical dependencies\n")
    print(vapply(c("rsvd", "Matrix"), function(p) as.character(packageVersion(p)), character(1)))
    cat("\nRNG\n")
    print(RNGkind())
    cat("\nThread environment\n")
    print(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                       "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")))
  }, file = environment_file)
}
cat(case, method, "repetition", repetition, "elapsed", elapsed,
    "seconds; peak RSS", round(peak_rss_kib/1024), "MiB; iterations", iterations,
    "converged", converged, "\n")
