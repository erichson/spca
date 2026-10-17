# Check the preserved source archive, including the complete PDF/HTML manuals.
# Usage: Rscript tools/release/check-source.R ARCHIVE NEW_OUTPUT_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply ARCHIVE and NEW_OUTPUT_DIRECTORY")
archive <- normalizePath(args[1L], mustWork = TRUE)
output <- args[2L]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                          no.. = TRUE))) {
  stop("Preserve existing diagnostics and choose an empty output directory")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
stopifnot(requireNamespace("rcmdcheck", quietly = TRUE),
          requireNamespace("V8", quietly = TRUE),
          requireNamespace("testthat", quietly = TRUE),
          nzchar(Sys.which("tidy")))
checksum <- function() {
  result <- system2("sha256sum", shQuote(archive), stdout = TRUE)
  if (!is.null(attr(result, "status"))) stop("SHA-256 calculation failed")
  result
}
before <- checksum()
writeLines(before, file.path(output, "archive-sha256.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
utils::write.table(installed.packages()[, c("Package", "Version", "Built")],
                   file.path(output, "package-versions.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
check_env <- c(NOT_CRAN = "false", `_R_CHECK_CRAN_INCOMING_` = "true",
               `_R_CHECK_CRAN_INCOMING_REMOTE_` = "true",
               `_R_CHECK_FORCE_SUGGESTS_` = "true",
               `_R_CHECK_SYSTEM_CLOCK_` = "true")
recorded_env <- c(check_env, Sys.getenv(c(
  "R_LIBS_USER", "PATH", "TEXMFHOME", "TEXMFCONFIG", "TEXMFVAR",
  "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS")))
writeLines(capture.output(print(recorded_env)),
           file.path(output, "check-environment.txt"))
writeLines(c("R CMD check --as-cran", archive), file.path(output, "command.txt"))
result <- rcmdcheck::rcmdcheck(archive, args = "--as-cran", check_dir = output,
                              env = check_env, error_on = "never", timeout = 1500)
saveRDS(result, file.path(output, "check-result.rds"))
writeLines(capture.output(print(result)), file.path(output, "check-summary.txt"))
stopifnot(identical(before, checksum()))
cat("ERRORS", length(result$errors), "WARNINGS", length(result$warnings),
    "NOTES", length(result$notes), "\n")
if (length(result$errors) || length(result$warnings)) quit(status = 1L)
