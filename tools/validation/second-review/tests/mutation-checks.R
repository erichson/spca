library(sparsepca)
library(testthat)
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[1] else ".")
output <- if (length(args) > 1L) args[2] else tempfile("sparsepca-mutations-")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
test_file_path <- file.path(repo, "tests/testthat/test-numerical-safeguards.R")
ns <- asNamespace("sparsepca")
loss <- get(".spca_loss", envir = ns)
progress <- get(".spca_progress", envir = ns)
result <- get(".spca_result", envir = ns)
change_body <- function(fun, old, new) {
  code <- paste(deparse(body(fun)), collapse = "\n")
  stopifnot(grepl(old, code, fixed = TRUE))
  body(fun) <- parse(text = gsub(old, new, code, fixed = TRUE))[[1]]
  fun
}
mutations <- list(
  no_loss_fallback = list(name = ".spca_loss", original = loss,
    changed = function(C, CA, CB, A, energy) {
      0.5 * (energy - 2 * sum(CA * CB) + sum(CB^2))
    }),
  increased_objective_converged = list(name = ".spca_progress", original = progress,
    changed = change_body(progress, "converged = FALSE, stop = TRUE",
                           "converged = TRUE, stop = TRUE")),
  ignore_evaluation_scale = list(name = ".spca_progress", original = progress,
    changed = change_body(progress, "abs(evaluation_scale)", "0")),
  forget_variance_pivot = list(name = ".spca_result", original = result,
    changed = change_body(result, "adjusted[qr_scores$pivot[keep]]",
                           "adjusted[keep]"))
)
run_one <- function(label, mutation = NULL) {
  if (!is.null(mutation)) {
    assignInNamespace(mutation$name, mutation$changed, ns = "sparsepca")
    on.exit(assignInNamespace(mutation$name, mutation$original, ns = "sparsepca"))
  }
  tests <- test_file(test_file_path, reporter = "silent")
  saveRDS(tests, file.path(output, paste0(label, "-results.rds")))
  tab <- as.data.frame(tests)
  tab$result <- NULL
  write.csv(tab, file.path(output, paste0(label, "-results.csv")), row.names = FALSE)
  data.frame(case = label, passed = sum(tab$passed), failed = sum(tab$failed),
              errors = sum(tab$error), warnings = sum(tab$warning),
              skipped = sum(tab$skipped))
}
summary <- run_one("unmodified")
for (label in names(mutations)) summary <- rbind(summary, run_one(label, mutations[[label]]))
stopifnot(summary$failed[1] == 0, summary$errors[1] == 0,
          all(summary$failed[-1] > 0), all(summary$errors == 0))
write.csv(summary, file.path(output, "mutation-summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
print(summary)
