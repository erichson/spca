#!/usr/bin/env Rscript
# Run generated help examples and executable README examples from a checkout.
# Installation instructions in README block 1 are intentionally not executed.
# Usage: Rscript tools/documentation/check-examples.R [root] [output_directory]
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1L]) else normalizePath(".")
out <- if (length(args) > 1L) args[2L] else
  file.path(root, "tools", "documentation", "results", "example-check")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
pkgload::load_all(root, quiet = TRUE)
results <- list()
run_code <- function(code, label, environment) {
  warnings <- character()
  failure <- NULL
  log <- capture.output({
    tryCatch(withCallingHandlers({
      for (expression in parse(text = code)) {
        visible <- withVisible(eval(expression, envir = environment))
        if (visible$visible) print(visible$value)
      }
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) failure <<- conditionMessage(e))
  })
  writeLines(c(log, if (length(warnings)) paste("WARNING:", warnings),
                 if (!is.null(failure)) paste("ERROR:", failure)),
               file.path(out, paste0(label, ".log")))
  results[[length(results) + 1L]] <<- data.frame(
    example = label,
    status = if (!is.null(failure)) "error" else if (length(warnings)) "warning" else "pass",
    warnings = paste(warnings, collapse = " | "),
    error = if (is.null(failure)) "" else failure,
    stringsAsFactors = FALSE
  )
}

rd_files <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
# Build local links from the source aliases. Installation normally supplies
# this map and evaluates build/install-time Rd macros such as DOI links.
html_links <- character()
for (rd in rd_files) {
  parsed <- tools::parse_Rd(rd)
  aliases <- vapply(Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"),
                           parsed), paste, collapse = "", FUN.VALUE = character(1))
  html_links[aliases] <- sub("[.]Rd$", ".html", basename(rd))
}
for (rd in rd_files) {
  label <- tools::file_path_sans_ext(basename(rd))
  code_path <- file.path(out, paste0(label, "-examples.R"))
  tools::Rd2ex(rd, out = code_path)
  if (file.exists(code_path)) {
    run_code(readLines(code_path, warn = FALSE), paste0("help-", label),
               new.env(parent = globalenv()))
  } else {
    results[[length(results) + 1L]] <- data.frame(
      example = paste0("help-", label), status = "no-examples",
      warnings = "", error = "", stringsAsFactors = FALSE)
  }
  tools::Rd2HTML(rd, out = file.path(out, paste0(label, ".html")),
                  package = "sparsepca", Links = html_links,
                  stages = c("build", "install", "render"))
}

readme <- readLines(file.path(root, "README.md"), warn = FALSE)
blocks <- list()
inside <- FALSE
for (line in readme) {
  if (identical(line, "```r")) {
    inside <- TRUE
    code <- character()
  } else if (inside && identical(line, "```")) {
    blocks[[length(blocks) + 1L]] <- code
    inside <- FALSE
  } else if (inside) code <- c(code, line)
}
readme_environment <- new.env(parent = globalenv())
grDevices::png(file.path(out, "README-plot-%02d.png"), width = 900, height = 600)
for (i in seq_along(blocks)) {
  if (any(grepl("install[.]packages|install_github", blocks[[i]]))) {
    results[[length(results) + 1L]] <- data.frame(
      example = paste0("README-", i), status = "skip-installation",
      warnings = "", error = "", stringsAsFactors = FALSE)
  } else run_code(blocks[[i]], paste0("README-", i), readme_environment)
}
grDevices::dev.off()
result <- do.call(rbind, results)
utils::write.csv(result, file.path(out, "example-results.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
print(result, row.names = FALSE)
if (any(result$status == "error")) stop("Documentation examples failed; inspect the saved logs.")
