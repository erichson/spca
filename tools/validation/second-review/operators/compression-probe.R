library(Matrix)
library(sparsepca)
for (f in list.files('/workspace/scratch/4668db2fbf57/spca/R', pattern = '[.]R$', full.names = TRUE)) source(f)
options(width = 160, warn = 1)
# Four disjoint pairs of positive/negative rows give a sparse zero-mean
# left factor; disjoint feature blocks yield exactly rank-four sparse input.
n <- 120L; p <- 72L; k <- 4L
left <- matrix(0, n, k)
right <- matrix(0, p, k)
for (j in seq_len(k)) {
  left[(2 * j - 1):(2 * j), j] <- c(1, -1) / sqrt(2)
  rows <- seq.int((j - 1) * 18 + 1, j * 18)
  right[rows, j] <- seq_along(rows) / sqrt(sum(seq_along(rows)^2))
}
exact <- left %*% diag(c(30, 15, 7, 3)) %*% t(right)
X <- Matrix(exact, sparse = TRUE)
cat('Exact case density:', nnzero(X) / prod(dim(X)), '\n')
for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
  target <- base::scale(exact, center = center, scale = scale)
  for (seed in c(8, 21, 79)) for (q in c(0, 2)) {
    set.seed(seed)
    fit <- rspca(X, k, alpha = 0, beta = 0, center = center, scale = scale,
                  o = 3, q = q, max_iter = 20, tol = 1e-12, verbose = FALSE)
    cat('exact center', center, 'scale', scale, 'seed', seed, 'q', q,
        'relative residual', norm(target - tcrossprod(fit$scores, fit$transform), 'F') / norm(target, 'F'),
        'score centered err', if(center) max(abs(colMeans(fit$scores))) else NA, '\n')
  }
}
# Sparse diagonal data supply a full-rank geometric spectrum. Embed enough
# zero rows that centering is a real operator correction, not a no-op.
spectrum <- exp(seq(log(30), log(0.03), length.out = p))
y <- matrix(0, n, p)
y[cbind(seq_len(p), seq_len(p))] <- spectrum
Y <- Matrix(y, sparse = TRUE)
for (center in c(FALSE, TRUE)) {
  target <- base::scale(y, center = center, scale = FALSE)
  d <- svd(target, nu = 0, nv = 0)$d
  tail2 <- sum(d[-seq_len(k)]^2)
  for (seed in c(8, 21, 79, 235, 417)) for (q in c(0, 1, 2, 4)) {
    set.seed(seed)
    fit <- rspca(Y, k, alpha = 0, beta = 0, center = center, o = 3, q = q,
                  max_iter = 20, tol = 1e-12, verbose = FALSE)
    error2 <- sum((target - tcrossprod(fit$scores, fit$transform))^2)
    cat('tail center', center, 'seed', seed, 'q', q, 'ratio', error2 / tail2,
        'approx model obj ratio', sum(fit$objective) / error2, '\n')
  }
}
cat('\nCompressed sparse/dense same-seed reconstruction quality comparison:\n')
for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
  target <- base::scale(y, center = center, scale = scale)
  optimal <- sum(svd(target, nu = 0, nv = 0)$d[-seq_len(k)]^2)
  fits <- lapply(list(y, Y), function(a) {set.seed(235); rspca(a, k, alpha = 0, beta = 0, center = center, scale = scale, o = 3, q = 2, max_iter = 20, tol = 1e-12, verbose = FALSE)})
  errors <- vapply(fits, function(f) sum((target - tcrossprod(f$scores, f$transform))^2) / optimal, 1)
  cat('center', center, 'scale', scale, 'ratios', paste(signif(errors, 12), collapse = '/'), '\n')
}
