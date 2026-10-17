library(Matrix)
library(sparsepca)
repo <- '/workspace/scratch/4668db2fbf57/spca'
for (f in list.files(file.path(repo, 'R'), pattern = '[.]R$', full.names = TRUE)) source(f)
options(width = 130, warn = 1)
cat('Runtime:', R.version.string, '; Matrix', as.character(packageVersion('Matrix')), '\n')
check_products <- function(sparse, label) {
  cat('\nCASE:', label, 'class:', class(sparse), '\n')
  original <- sparse
  x <- as.matrix(sparse)
  missing_rows <- rowSums(is.na(x)) > 0
  x <- x[!missing_rows, , drop = FALSE]
  for (center in c(FALSE, TRUE)) for (scale in c(FALSE, TRUE)) {
    prep <- suppressWarnings(.spca_preprocess(sparse, center, scale, sparse = TRUE))
    target <- x
    expected_center <- if (center) colMeans(target) else FALSE
    if (center) target <- sweep(target, 2, expected_center, '-')
    expected_scale <- if (scale) apply(target, 2, function(v) sqrt(sum(v * v)/(nrow(x)-1))) else FALSE
    if (scale) {
      expected_scale[expected_scale < 1e-8] <- 1
      target <- sweep(target, 2, expected_scale, '/')
    }
    right <- matrix(sin(seq_len(ncol(x) * 3)), ncol(x), 3)
    left <- matrix(cos(seq_len(nrow(x) * 3)), nrow(x), 3)
    errors <- c(forward = max(abs(prep$forward(right) - target %*% right)),
                adjoint = max(abs(prep$adjoint(left) - crossprod(target, left))),
                energy = abs(prep$energy - sum(target * target)))
    cat('center', center, 'scale', scale, ':', paste(names(errors), signif(errors, 6), collapse = '; '), '\n')
    stopifnot(isTRUE(all.equal(prep$center, expected_center, tolerance = 1e-12)),
              isTRUE(all.equal(prep$scale, expected_scale, tolerance = 1e-12)),
              max(abs(prep$forward(right) - target %*% right)) <= 1e-10 * max(1, abs(target %*% right)),
              max(abs(prep$adjoint(left) - crossprod(target, left))) <= 1e-10 * max(1, abs(crossprod(target, left))),
              abs(prep$energy - sum(target * target)) <= 1e-10 * max(1, sum(target * target)),
              identical(original, sparse))
  }
}
set.seed(771)
x <- matrix(rnorm(23 * 7), 23, 7)
x[abs(x) < 1.1] <- 0
x[, 5] <- 0
x[, 6] <- 3
x[, 7] <- seq_len(23) * 1e-10
rownames(x) <- paste0('obs', seq_len(nrow(x)))
colnames(x) <- paste0('feature', seq_len(ncol(x)))
check_products(Matrix(x, sparse = TRUE), 'mixed sparse/constant/tiny columns')
explicit_zero <- as(Matrix(matrix(0, 23, 7), sparse = TRUE), 'generalMatrix')
explicit_zero <- sparseMatrix(i = rep(seq_len(23), 7), j = rep(seq_len(7), each = 23), x = as.vector(x), dims = dim(x), dimnames = dimnames(x))
check_products(explicit_zero, 'explicit stored zeros')
sym <- forceSymmetric(rsparsematrix(8, 8, density = 0.3))
sym[2, 5] <- NA
sym[7, 7] <- NA
check_products(sym, 'symmetric missing offdiagonal and diagonal')
tri <- as(Diagonal(8), 'triangularMatrix')
tri[2, 5] <- NA
check_products(tri, 'unit triangular missing stored entry')
dup <- sparseMatrix(i = c(1,1,2,2,3,4,6), j = c(1,1,2,2,3,1,3), x = c(2,-2,1,NA,3,4,-1), dims = c(9,4), repr = 'T')
check_products(dup, 'triplet duplicates plus missing row')
large <- matrix(1e150, 20, 5)
check_products(Matrix(large, sparse = TRUE), 'constant large offset')
cat('\nPredictions with zero rows:\n')
fit <- spca(x, 2, scale = TRUE, verbose = FALSE)
for (new in list(x[FALSE,,drop = FALSE], Matrix(x[FALSE,,drop = FALSE], sparse = TRUE), as(Matrix(x[FALSE,,drop = FALSE], sparse = TRUE), 'TsparseMatrix'), as(Matrix(x[FALSE,,drop = FALSE], sparse = TRUE), 'RsparseMatrix'))) {
  pred <- predict(fit, new)
  cat(class(new)[1], '->', paste(dim(pred), collapse = 'x'), '; names', paste(names(dimnames(pred)), collapse = '/'), '\n')
  stopifnot(identical(dim(pred), c(0L, 2L)))
}
cat('\nIncomplete/duplicate training and prediction names:\n')
for (names in list(c(NA, paste0('f', 2:7)), c('', paste0('f', 2:7)), rep('f', 7))) {
  named <- x
  colnames(named) <- names
  fit <- spca(named, 2, scale = TRUE, verbose = FALSE)
  value <- tryCatch(predict(fit, named[, 7:1, drop = FALSE]), error = function(e) e)
  if (inherits(value, 'error')) cat('names:', paste(names, collapse = ','), '-> ERROR:', conditionMessage(value), '\n') else cat('names:', paste(names, collapse = ','), '-> accepted, scores max error', max(abs(value - fit$scores)), '\n')
}
cat('\nExtreme scaling behavior (not assumed a correctness bug):\n')
for (multiplier in c(1e-160, 1e-150, 1e140, 1e150, 1e153, 1e154)) {
  y <- matrix(c(-1, 0, 1, 2, 1, -3), 3, 2) * multiplier
  for (is_sparse in c(FALSE, TRUE)) {
    value <- tryCatch(.spca_preprocess(if (is_sparse) Matrix(y, sparse = TRUE) else y, TRUE, TRUE, sparse = is_sparse), error = function(e) e)
    cat(multiplier, if(is_sparse) 'sparse' else 'dense', ':', if (inherits(value, 'error')) conditionMessage(value) else paste('energy', value$energy), '\n')
  }
}
cat('\nALL explicit operator comparisons passed.\n')
true_symmetric <- sparseMatrix(i = c(1, 2, 4, 3), j = c(1, 5, 8, 3),
                              x = c(2, NA, 4, 3), dims = c(8, 8), symmetric = TRUE)
stopifnot(inherits(true_symmetric, 'symmetricMatrix'))
check_products(true_symmetric, 'structurally symmetric missing offdiagonal')
