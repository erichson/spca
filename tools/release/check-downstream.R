# Run current downstream source archives against a selected sparsepca library.
# This development script is excluded from the CRAN source package.
#
# Usage (run once with baseline and once with candidate R_LIBS_USER):
# Rscript tools/release/check-downstream.R ARCHIVE_DIR OUTPUT_DIR VERSION \
#   categoryEncodings clusterability robCompositions SparseBiplots
#
# Requires rcmdcheck. The supplied archives, runtime, BLAS/thread settings,
# optional dependencies, and check flags should be identical for both runs.
# Unlike a CRAN submission preflight, these checks do not run incoming checks:
# the already-published downstream packages are not being submitted again.
# PDF manuals are excluded; installation, code analysis, examples, and tests run.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: check-downstream.R ARCHIVE_DIR OUTPUT_DIR SPARSEPCA_VERSION [PACKAGES]")
}
archive_dir <- normalizePath(args[1L], mustWork = TRUE)
output_dir <- args[2L]
expected_version <- args[3L]
packages <- if (length(args) > 3L) args[-seq_len(3L)] else {
  c("categoryEncodings", "clusterability", "robCompositions", "SparseBiplots")
}
stopifnot(requireNamespace("rcmdcheck", quietly = TRUE),
          as.character(utils::packageVersion("sparsepca")) == expected_version)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

flags <- "--no-manual"
check_env <- c(NOT_CRAN = "false", `_R_CHECK_CRAN_INCOMING_` = "false",
               `_R_CHECK_FORCE_SUGGESTS_` = "true",
               OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
write.table(utils::installed.packages()[, c("Package", "Version", "Built")],
            file.path(output_dir, "package-versions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("sparsepca:", expected_version),
             paste("library:", find.package("sparsepca")),
             paste("flags:", paste(flags, collapse = " ")),
             paste(names(check_env), check_env, sep = "=")),
           file.path(output_dir, "check-settings.txt"))

summaries <- vector("list", length(packages))
names(summaries) <- packages
for (package in packages) {
  archive <- list.files(archive_dir,
    pattern = paste0("^", package, "_[0-9.-]+[.]tar[.]gz$"), full.names = TRUE)
  if (length(archive) != 1L) stop("Expected one source archive for ", package)
  cat("\nChecking", basename(archive), "with sparsepca", expected_version, "\n")
  package_dir <- file.path(output_dir, package)
  dir.create(package_dir, showWarnings = FALSE)
  write.table(data.frame(archive = basename(archive),
                         md5 = unname(tools::md5sum(archive))),
              file.path(package_dir, "archive-identity.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  result <- rcmdcheck::rcmdcheck(archive, args = flags, check_dir = package_dir,
                               error_on = "never", timeout = 1800,
                               env = check_env)
  saveRDS(result, file.path(package_dir, "check-result.rds"))
  writeLines(capture.output(print(result)), file.path(package_dir, "check-summary.txt"))
  summaries[[package]] <- data.frame(
    package = package, version = result$version,
    sparsepca = expected_version, errors = length(result$errors),
    warnings = length(result$warnings), notes = length(result$notes),
    stringsAsFactors = FALSE)
  utils::write.csv(do.call(rbind, summaries[!vapply(summaries, is.null, logical(1))]),
                   file.path(output_dir, "summary.csv"), row.names = FALSE)
}

