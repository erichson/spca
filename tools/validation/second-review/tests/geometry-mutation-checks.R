library(sparsepca)
library(testthat)
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[1] else ".")
output <- if (length(args) > 1L) args[2] else tempfile("sparsepca-mutations-")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source_path <- file.path(repo, "tests/testthat/test-objective-geometry.R")
cases <- parse(source_path, keep.source = FALSE)
# Limit fault-injection to the two mathematical oracle blocks under review.
focused_path <- file.path(output, "focused-objective-geometry.R")
writeLines(vapply(cases[1:2], function(x) paste(deparse(x), collapse = "\n"),
                   character(1)), focused_path)
ns <- asNamespace("sparsepca")
solver <- get(".spca_solve", envir = ns)
change_body <- function(fun, old, new) {
  code <- paste(deparse(body(fun)), collapse = "\n")
  stopifnot(grepl(old, code, fixed = TRUE))
  body(fun) <- parse(text = gsub(old, new, code, fixed = TRUE))[[1]]
  fun
}
mutations <- list(
  omit_ridge_shrinkage = change_body(solver, "shrink <- 1 + controls$beta", "shrink <- 1"),
  reverse_smooth_gradient = change_body(solver, "CA - CB", "CB - CA"),
  treat_alpha_as_absolute = change_body(solver, "threshold <- controls$alpha/shrink",
    "threshold <- controls$alpha/(shrink * factor$norm2)")
)
run_one <- function(label, changed = NULL) {
  if (!is.null(changed)) {
    assignInNamespace(".spca_solve", changed, ns = "sparsepca")
    on.exit(assignInNamespace(".spca_solve", solver, ns = "sparsepca"))
  }
  tests <- test_file(focused_path, reporter = "silent")
  saveRDS(tests, file.path(output, paste0("geometry-", label, "-results.rds")))
  tab <- as.data.frame(tests)
  tab$result <- NULL
  write.csv(tab, file.path(output, paste0("geometry-", label, "-results.csv")),
              row.names = FALSE)
  data.frame(case = label, block = seq_len(nrow(tab)), passed = tab$passed,
              failed = tab$failed, errors = tab$error, warnings = tab$warning,
              skipped = tab$skipped)
}
summary <- run_one("unmodified")
for (label in names(mutations)) summary <- rbind(summary, run_one(label, mutations[[label]]))
write.csv(summary, file.path(output, "geometry-mutation-summary.csv"), row.names = FALSE)
print(summary)
stopifnot(all(summary$failed[summary$case == "unmodified"] == 0),
          all(summary$errors == 0),
          all(summary$failed[summary$case != "unmodified"] > 0))

# Validate the oracle's smoothness and distance from a hard-threshold tie.
set.seed(802)
X <- matrix(rnorm(18 * 6), 18, 6)
X[, 2:6] <- X[, 2:6] + 0.4 * X[, 1:5]
B <- matrix(rnorm(12), 6, 2) / 3
value <- function(B) {
  XB <- X %*% B
  0.5 * sum(X^2) + 0.5 * sum(XB^2) -
      sum(svd(crossprod(X, XB), nu = 0, nv = 0)$d)
}
h <- 1e-5
gradient <- B
for (j in seq_along(B)) {
  plus <- minus <- B
  plus[j] <- plus[j] + h
  minus[j] <- minus[j] - h
  gradient[j] <- (value(plus) - value(minus))/(2*h)
}
spectral2 <- svd(X, nu = 0, nv = 0)$d[1]^2
alpha <- 0.035
beta <- 0.2
step <- 1 / (spectral2 * (1 + beta))
z <- B - step * (gradient + beta * spectral2 * B)
metadata <- c(min_cross_singular_value = min(svd(crossprod(X, X %*% B),
                                                nu=0,nv=0)$d),
               min_hard_cost_margin = min(abs(0.5*z^2 - step*alpha*spectral2)),
               largest_data_singular_value_squared = spectral2)
for (penalty in c("elasticnet", "l0")) {
  expected <- if (penalty == "elasticnet") {
    sign(z) * pmax(abs(z) - step*alpha*spectral2, 0)
  } else z * (0.5*z^2 > step*alpha*spectral2)
  fit <- spca(X, 2, alpha=alpha,beta=beta,center=FALSE,penalty=penalty,
               B_init=B,max_iter=1,verbose=FALSE)
  metadata[paste0("max_loading_error_", penalty)] <- max(abs(fit$loadings - expected))
}
write.csv(data.frame(metric=names(metadata), value=unname(metadata)),
          file.path(output, "geometry-oracle-diagnostics.csv"), row.names=FALSE)
print(metadata)
