target_test_matrix <- function() {
  rbind(diag(c(4, 3, 2)), -diag(c(4, 3, 2)))
}

# Base tracing keeps this instrumentation compatible with testthat 3.0.0.
# Always restore the namespace, including when the evaluated call fails.
target_count_calls <- function(code) {
  calls <- new.env(parent = emptyenv())
  functions <- c(".spca_preprocess", ".spca_factor", ".spca_result")
  namespace <- asNamespace("sparsepca")
  traced <- character()
  on.exit(for (name in traced) {
    suppressMessages(untrace(name, where = namespace))
  }, add = TRUE)
  for (name in functions) {
    calls[[name]] <- 0L
    tracer <- substitute(assign(KEY, get(KEY, envir = COUNTER) + 1L,
                                envir = COUNTER),
                          list(COUNTER = calls, KEY = name))
    suppressMessages(trace(name, tracer = tracer, print = FALSE,
                           where = namespace))
    traced <- c(traced, name)
  }
  value <- force(code)
  list(value = value, calls = vapply(functions, function(name) calls[[name]],
                                    integer(1)))
}

test_that("selection handles nonmonotone counts and only eligible fits", {
  pick <- sparsepca:::.spca_target_pick
  expect_identical(pick(c(12L, 8L, 10L, 7L, 9L),
                         c(TRUE, TRUE, FALSE, TRUE, TRUE), 9L), 5L)
  # An unconverged exact count cannot defeat a converged approximate count.
  expect_identical(pick(c(2L, 5L, 3L, 4L),
                         c(TRUE, TRUE, TRUE, FALSE), 4L), 3L)
  # Prefer fewer variables on an equal-distance tie, then the first fit.
  expect_identical(pick(c(7L, 11L, 7L), rep(TRUE, 3), 9L), 1L)
  expect_identical(pick(c(5L, 3L, 3L), rep(TRUE, 3), 4L), 2L)
  expect_identical(pick(c(3L, 4L), c(FALSE, FALSE), 4L), NA_integer_)
})

test_that("refinement uses observed crossings in either direction", {
  next_alpha <- sparsepca:::.spca_target_next
  alpha <- c(0, 1, 2, 100)
  # The widest interval has no crossing. Both rising and falling crossings
  # take priority even though all intervals touch a count equally near 5.
  for (counts in list(c(4, 6, 4, 4), c(6, 4, 6, 6))) {
    value <- next_alpha(alpha, counts, 5, 1)
    expect_true(value > 0 && value < 2)
    expect_false(value %in% alpha)
  }
  expect_equal(next_alpha(c(1, 9), c(4, 6), 5, 1), 3)
  expect_equal(next_alpha(c(0, 2), c(6, 4), 5, 1), 1)
  # A plateau still permits exploration; the target need not be bracketed.
  value <- next_alpha(c(0, 1, 4), c(6, 6, 6), 5, 1)
  expect_true(is.finite(value) && value > 0 && value < 4)
  expect_false(value %in% c(0, 1, 4))
})

test_that("refinement handles repeated and extreme penalties safely", {
  next_alpha <- sparsepca:::.spca_target_next
  expect_identical(next_alpha(c(.1, .1), c(2, 4), 3, .1), NA_real_)
  expect_equal(next_alpha(c(0, 2, 2), c(8, 4, 6), 5, 2), 1)
  # There is no representable floating-point value strictly between these.
  expect_identical(next_alpha(c(1, 1 + .Machine$double.eps),
                                c(2, 4), 3, 1), NA_real_)
  alpha <- c(0, 1e-300, 1e300)
  value <- next_alpha(alpha, c(9, 6, 2), 5, 1e-300)
  expect_true(is.finite(value) && value > 1e-300 && value < 1e300)
  expect_false(value %in% alpha)
})

test_that("a feature count is the union of nonzero loading rows", {
  X <- outer(-5:5, c(1, 2, 3, 0))
  initial <- cbind(c(1, 2, 3, 0), c(1, 2, 3, 0)) / sqrt(14)
  colnames(X) <- paste0("feature", seq_len(ncol(X)))
  answer <- spca_target(X, 2, 3, alpha = 0, beta = .1, B_init = initial,
                        tol = 1e-10, max_evals = 3, verbose = FALSE)
  expect_s3_class(answer, "spca_target")
  expect_s3_class(answer$fit, "spca")
  B <- answer$fit$loadings
  expect_identical(answer$selected_features, 3L)
  expect_equal(answer$selected_features, sum(rowSums(B != 0) > 0))
  expect_gt(sum(B != 0), answer$selected_features)
  expect_equal(unname(answer$nonzero[1, ]), unname(colSums(B != 0)))
  expect_identical(dim(answer$nonzero), c(1L, 2L))
  expect_identical(colnames(answer$nonzero), c("PC1", "PC2"))
  expect_identical(answer$fit$feature_names, colnames(X))
  expect_equal(predict(answer$fit, X), answer$fit$scores)
  expect_identical(dim(summary(answer$fit)), c(4L, 2L))
})

test_that("supplied penalties keep their order and stop at the first match", {
  X <- target_test_matrix()
  for (penalty in c("elasticnet", "l0")) {
    answer <- spca_target(X, 1, 0, alpha = c(.05, 1, .2), beta = .1,
                          penalty = penalty, tol = 1e-10, verbose = FALSE)
    expect_identical(answer$search$alpha, c(.05, 1))
    expect_identical(answer$search$warm_start, 0:1)
    expect_identical(answer$search$selected_features, c(1L, 0L))
    expect_true(answer$target_met)
    expect_identical(answer$index, 2L)
    expect_identical(answer$alpha, 1)
    expect_identical(answer$stop_reason, "target_met")
    expect_true(all(answer$fit$loadings == 0))
    expect_identical(answer$penalty, penalty)
  }
  duplicate <- spca_target(X, 1, 0, alpha = c(.05, .05, 1), beta = .1,
                           tol = 1e-10, verbose = FALSE)
  expect_identical(duplicate$search$alpha, c(.05, .05, 1))
  expect_identical(duplicate$search$warm_start, 0:2)
})

test_that("count tolerance accepts a nearby count explicitly", {
  answer <- spca_target(target_test_matrix(), 1, 2, alpha = c(.05, 1),
                        count_tol = 1, beta = .1, verbose = FALSE)
  expect_true(answer$target_met)
  expect_identical(answer$n_features, 2L)
  expect_identical(answer$selected_features, 1L)
  expect_identical(answer$count_tol, 1L)
  expect_identical(answer$search$difference, -1L)
  expect_identical(nrow(answer$search), 1L)
  expect_identical(answer$stop_reason, "target_met")
})

test_that("missed counts return the closest converged fit within the budget", {
  # A diagonal covariance and k = 1 keep this path on one coordinate. The
  # requested two variables are skipped even though they are a valid target.
  for (penalty in c("elasticnet", "l0")) {
    answer <- spca_target(target_test_matrix(), 1, 2,
                          alpha = c(.05, 1, .2), beta = .1,
                          penalty = penalty, max_evals = 6,
                          tol = 1e-10, verbose = FALSE)
    expect_false(answer$target_met)
    expect_identical(answer$selected_features, 1L)
    expect_identical(answer$index, 1L)
    expect_identical(nrow(answer$search), 6L)
    expect_identical(answer$stop_reason, "max_evals")
    expect_true(all(answer$search$converged))
    expect_identical(answer$search$alpha[1:3], c(.05, 1, .2))
    expect_identical(length(unique(answer$search$alpha)), 6L)
    expect_true(all(answer$search$warm_start >= 0L &
                      answer$search$warm_start < seq_len(nrow(answer$search))))
  }
})

test_that("automatic search honors small budgets and cannot refine one value", {
  X <- target_test_matrix()
  one <- spca_target(X, 1, 2, max_evals = 1, verbose = FALSE)
  expect_identical(one$search$alpha, 0)
  expect_identical(one$stop_reason, "max_evals")
  two <- spca_target(X, 1, 2, max_evals = 2, verbose = FALSE)
  expect_identical(two$search$alpha, c(0, 1))
  expect_identical(two$stop_reason, "max_evals")
  single <- spca_target(X, 1, 2, alpha = .1, max_evals = 4, verbose = FALSE)
  expect_identical(nrow(single$search), 1L)
  expect_identical(single$stop_reason, "no_refinement")
  expect_false(single$target_met)
})

test_that("an unconverged exact count does not produce a selected fit", {
  expect_warning(answer <- spca_target(target_test_matrix(), 1, 0,
                                       alpha = c(1, 2), max_evals = 2,
                                       max_iter = 1, verbose = FALSE),
                  "No candidate fit converged")
  expect_identical(answer$search$selected_features, c(0L, 0L))
  expect_false(any(answer$search$converged))
  expect_null(answer$fit)
  expect_identical(answer$index, NA_integer_)
  expect_identical(answer$selected_features, NA_integer_)
  expect_identical(answer$alpha, NA_real_)
  expect_false(answer$target_met)
  expect_identical(answer$stop_reason, "no_converged_fit")
  output <- capture.output(printed <- withVisible(print(answer)))
  expect_match(paste(output, collapse = "\n"), "No converged candidate")
  expect_false(printed$visible)
  expect_identical(printed$value, answer)
})

test_that("zero energy and deficient rank retain valid fit dimensions", {
  zero <- spca_target(matrix(0, 6, 4), 2, 0, max_iter = 1,
                      max_evals = 3, verbose = FALSE)
  expect_true(zero$target_met)
  expect_identical(zero$selected_features, 0L)
  expect_true(zero$fit$converged)
  expect_equal(zero$fit$objective, 0)
  expect_identical(dim(zero$fit$loadings), c(4L, 2L))
  expect_identical(dim(zero$fit$scores), c(6L, 2L))
  expect_true(all(zero$fit$scores == 0))
  missed <- spca_target(matrix(0, 6, 4), 2, 1, max_evals = 3, verbose = FALSE)
  expect_false(missed$target_met)
  expect_identical(missed$index, 1L)
  expect_identical(missed$stop_reason, "max_evals")
  X <- outer(seq_len(12), c(1, 1, 0, 0))
  deficient <- spca_target(X, 3, 2, alpha = .2, beta = .1,
                           tol = 1e-7, verbose = FALSE)
  expect_true(deficient$target_met)
  expect_identical(dim(deficient$fit$loadings), c(4L, 3L))
  expect_equal(crossprod(deficient$fit$transform), diag(3), tolerance = 1e-12)
  expect_true(all(is.finite(deficient$fit$objective)))
})

test_that("the selected fit replays its recorded warm-start ancestry", {
  X <- target_test_matrix()
  initial <- matrix(c(.7, 0, 0), 3, 1)
  for (method in c("spca", "rspca")) {
    set.seed(331)
    answer <- spca_target(X, 1, 2, alpha = c(.05, 1, .2), beta = .1,
                          method = method, o = 1, q = 1, B_init = initial,
                          max_evals = 6, tol = 1e-10, verbose = FALSE)
    replay <- vector("list", nrow(answer$search))
    fun <- if (method == "spca") spca else rspca
    for (i in seq_along(replay)) {
      ancestor <- answer$search$warm_start[i]
      B <- if (ancestor == 0L) initial else replay[[ancestor]]$loadings
      arguments <- list(X = X, k = 1, alpha = answer$search$alpha[i],
                        beta = .1, B_init = B, tol = 1e-10, verbose = FALSE)
      if (method == "rspca") arguments <- c(arguments, list(o = 1, q = 1))
      # Recreate the same sketch for every independent randomized reference.
      set.seed(331)
      replay[[i]] <- do.call(fun, arguments)
      expect_equal(sum(rowSums(replay[[i]]$loadings != 0) > 0),
                    answer$search$selected_features[i])
      expect_identical(replay[[i]]$converged, answer$search$converged[i])
      expect_identical(replay[[i]]$iterations, answer$search$iterations[i])
    }
    expect_identical(answer$search$warm_start[4:6], rep(1L, 3))
    selected <- replay[[answer$index]]
    expect_equal(answer$fit, selected, tolerance = 1e-12)
    expect_identical(answer$fit$alpha, selected$alpha)
  }
})

test_that("the search preprocesses once and constructs only one full result", {
  measured <- target_count_calls(spca_target(target_test_matrix(), 1, 2,
                                             alpha = c(.05, 1, .2),
                                             max_evals = 6, verbose = FALSE))
  expect_identical(nrow(measured$value$search), 6L)
  expect_identical(unname(measured$calls), rep(1L, 3))
})

test_that("refinements retain valid ancestry when no candidate converges", {
  X <- target_test_matrix()
  expect_warning(answer <- spca_target(X, 1, 2, alpha = c(.05, 1),
                                       beta = .1, max_iter = 1,
                                       max_evals = 4, verbose = FALSE),
                  "No candidate fit converged")
  expect_identical(answer$search$warm_start, c(0L, 1L, 1L, 1L))
  expect_false(any(answer$search$converged))
  expect_null(answer$fit)
  expect_identical(answer$stop_reason, "no_converged_fit")
  replay <- vector("list", nrow(answer$search))
  for (i in seq_along(replay)) {
    ancestor <- answer$search$warm_start[i]
    initial <- if (ancestor == 0L) NULL else replay[[ancestor]]$loadings
    replay[[i]] <- spca(X, 1, alpha = answer$search$alpha[i], beta = .1,
                        B_init = initial, max_iter = 1, verbose = FALSE)
    expect_equal(sum(rowSums(replay[[i]]$loadings != 0) > 0),
                  answer$search$selected_features[i])
    expect_equal(unname(colSums(replay[[i]]$loadings != 0)),
                  unname(answer$nonzero[i, ]))
  }
})

test_that("randomized searches consume one sketch and respect the caller RNG", {
  set.seed(9)
  X <- matrix(rnorm(40 * 15), 40, 15)
  set.seed(910)
  before <- .Random.seed
  rspca(X, 2, o = 3, q = 1, max_iter = 3, tol = 1, verbose = FALSE)
  after_one <- .Random.seed
  expect_false(identical(before, after_one))
  set.seed(910)
  answer <- spca_target(X, 2, 0, method = "rspca", o = 3, q = 1,
                        alpha = c(.001, .01, .05), max_evals = 3,
                        max_iter = 3, tol = 1, verbose = FALSE)
  expect_identical(nrow(answer$search), 3L)
  expect_identical(.Random.seed, after_one)
  seed <- .Random.seed
  spca_target(X, 2, 0, alpha = .01, max_evals = 1,
              max_iter = 3, tol = 1, verbose = FALSE)
  expect_identical(.Random.seed, seed)
})

test_that("sparse and dense inputs give equivalent randomized selections", {
  X <- Matrix::bdiag(outer(seq_len(9), c(1, 2, 3)),
                     outer(seq_len(7), c(2, 5)),
                     outer(seq_len(5), c(1, 4, 7, 9)))
  run <- function(X) {
    set.seed(77)
    spca_target(X, 2, 7, method = "rspca", o = 2, q = 1,
                 alpha = c(.02, .1, .3), max_evals = 3,
                 max_iter = 1000, tol = 1e-7, verbose = FALSE)
  }
  sparse <- run(X)
  dense <- run(as.matrix(X))
  expect_s3_class(sparse$fit, "rspca")
  expect_identical(sparse$index, dense$index)
  expect_identical(sparse$search$selected_features, dense$search$selected_features)
  expect_identical(sparse$search$converged, dense$search$converged)
  expect_identical(sparse$nonzero, dense$nonzero)
  # Products are invariant to simultaneous signs or basis changes in factors.
  sparse_reconstruction <- tcrossprod(sparse$fit$scores, sparse$fit$transform)
  dense_reconstruction <- tcrossprod(dense$fit$scores, dense$fit$transform)
  expect_equal(sparse_reconstruction, dense_reconstruction, tolerance = 1e-6)
  expect_equal(sparse$fit$objective, dense$fit$objective, tolerance = 1e-6)
  expect_equal(predict(sparse$fit, X), sparse$fit$scores, tolerance = 1e-12)
})

test_that("feature targeting rejects invalid counts and search controls", {
  X <- target_test_matrix()
  for (value in list(-1, 1.5, NA_real_, Inf, 1 + 1i, "1", TRUE, c(1, 2), 4)) {
    expect_error(spca_target(X, 1, value, verbose = FALSE), "n_features")
  }
  for (value in list(-1, .5, NA_real_, Inf, .Machine$integer.max + 1)) {
    expect_error(spca_target(X, 1, 1, count_tol = value, verbose = FALSE),
                  "count_tol")
  }
  for (value in list(0, 1.5, NA_real_, Inf, .Machine$integer.max + 1)) {
    expect_error(spca_target(X, 1, 1, max_evals = value, verbose = FALSE),
                  "max_evals")
  }
  for (value in list(numeric(), -1, NA_real_, Inf, 1 + 1i, "0.1")) {
    expect_error(spca_target(X, 1, 1, alpha = value, verbose = FALSE), "alpha")
  }
  expect_error(spca_target(X, 1, 1, alpha = c(0, .1, 1), max_evals = 2),
                "max_evals")
  expect_error(spca_target(X, 1, 1, B_init = matrix(1, 2, 1)), "B_init")
})
