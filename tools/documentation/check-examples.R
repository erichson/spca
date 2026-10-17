#!/usr/bin/env Rscript
# Verify source help and labeled README examples against an installed package.
# Usage: Rscript --vanilla tools/documentation/check-examples.R ROOT OUTPUT [LIBRARY]
# No installation/download is performed. Run installation instructions manually.

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 3L) stop("Expected ROOT OUTPUT [LIBRARY].", call. = FALSE)
  root <- normalizePath(if (length(args)) args[1L] else ".", mustWork = TRUE)
  out <- if (length(args) > 1L) args[2L] else
    file.path(root, "tools", "documentation", "results", "example-check")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, mustWork = TRUE)
  if (length(list.files(out, all.files = TRUE, no.. = TRUE))) {
    stop("Use a fresh OUTPUT directory so old logs and plots cannot enter this run.",
         call. = FALSE)
  }
  for (subdir in c("code", "logs", "plots", "html")) {
    dir.create(file.path(out, subdir), showWarnings = FALSE)
  }
  if (length(args) == 3L) {
    requested <- strsplit(args[3L], .Platform$path.sep, fixed = TRUE)[[1L]]
    requested <- normalizePath(requested, mustWork = TRUE)
    .libPaths(unique(c(requested, .libPaths())))
  }
  description <- read.dcf(file.path(root, "DESCRIPTION"))[1L, ]
  package <- unname(description["Package"])
  package_path <- find.package(package, quiet = TRUE)
  if (!length(package_path)) {
    stop("Install ", package, " outside this runner, then set R_LIBS_USER or pass LIBRARY.",
         call. = FALSE)
  }
  installed_version <- as.character(utils::packageVersion(package))
  if (package_version(installed_version) != package_version(description["Version"])) {
    stop("The source documentation expects ", package, " ", description["Version"],
         "; the installed package is ", installed_version, " at ", package_path,
         ". Install the candidate separately or select its library.", call. = FALSE)
  }
  library(package, character.only = TRUE)
  # Exclude this runner's global helpers. Examples see installed exports and
  # ordinary attached R packages; no pkgload-only internals are exposed.
  example_parent <- as.environment(paste0("package:", package))
  previous_wd <- getwd()
  on.exit(setwd(previous_wd), add = TRUE)
  setwd(root)

  source_files <- sort(unique(c(
    file.path(root, c("DESCRIPTION", "NAMESPACE", "README.md")),
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
    list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE),
    list.files(file.path(root, "tools", "documentation"),
               pattern = "[.](R|py|md)$", full.names = TRUE)
  )))
  hashes <- function(paths) {
    data.frame(file = substring(paths, nchar(root) + 2L),
               md5 = unname(tools::md5sum(paths)), stringsAsFactors = FALSE)
  }
  source_before <- hashes(source_files)
  utils::write.csv(source_before, file.path(out, "source-manifest.csv"), row.names = FALSE)
  git_value <- function(arguments) {
    if (!nzchar(Sys.which("git"))) return(NA_character_)
    value <- suppressWarnings(system2("git", c("-C", shQuote(root), arguments),
                                       stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(value, "status"))) return(NA_character_)
    paste(value, collapse = "\n")
  }
  dependencies <- c(package, "Matrix", "rsvd")
  versions <- setNames(vapply(dependencies, function(name) {
    as.character(utils::packageVersion(name))
  }, character(1)), dependencies)
  manifest <- list(
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source_root = root, output_directory = out,
    source_version = unname(description["Version"]),
    git_commit = git_value(c("rev-parse", "HEAD")),
    git_tree = git_value(c("rev-parse", "HEAD^{tree}")),
    git_status = git_value(c("status", "--porcelain")),
    installed_package_path = package_path, installed_versions = versions,
    library_paths = .libPaths(), source_hash_algorithm = "MD5",
    flags = list(installed_package = TRUE, pkgload = FALSE,
                 warnings_fail = TRUE, help_environments = "independent",
                 README_environment = "sequential",
                 unlabeled_R_chunks_fail = TRUE,
                 unevaluated_chunks_require_explicit_markup = TRUE,
                 local_HTML_stylesheet_and_index = TRUE,
                 external_help_preview = "stats::prcomp",
                 external_help_examples_executed = FALSE),
    evidence_scope = paste("Source README/help examples executed against the identified installed package.",
                            "Matching version strings alone do not establish archive equivalence;",
                            "this check does not replace checking an exact release archive.")
  )
  installed_files <- c(file.path(package_path, "DESCRIPTION"),
                       list.files(file.path(package_path, "R"), full.names = TRUE))
  utils::write.csv(data.frame(file = installed_files,
                               md5 = unname(tools::md5sum(installed_files))),
                    file.path(out, "installed-package-manifest.csv"), row.names = FALSE)
  results <- list()
  record <- function(label, kind, status, start = NA_integer_, end = NA_integer_,
                     reason = "", warnings = character(), error = "") {
    results[[length(results) + 1L]] <<- data.frame(
      example = label, kind = kind, status = status,
      source_line_start = start, source_line_end = end, reason = reason,
      warnings = paste(warnings, collapse = " | "), error = error,
      stringsAsFactors = FALSE)
  }
  execute <- function(label, kind, action, start = NA_integer_, end = NA_integer_,
                       plot = FALSE) {
    warnings <- messages <- character()
    failure <- NULL
    device <- NULL
    log <- capture.output({
      tryCatch(withCallingHandlers({
        if (plot) {
          filename <- file.path(out, "plots", paste0(label, "-%03d"))
          if (capabilities("png")) {
            grDevices::png(paste0(filename, ".png"), width = 1200, height = 800,
                            res = 150)
          } else {
            grDevices::pdf(paste0(filename, ".pdf"), width = 8, height = 5.3,
                            onefile = FALSE)
          }
          device <- grDevices::dev.cur()
        }
        action()
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }, message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }), error = function(e) failure <<- conditionMessage(e), finally = {
        if (!is.null(device) && device %in% grDevices::dev.list()) {
          grDevices::dev.off(device)
        }
      })
    })
    writeLines(c(log, if (length(messages)) paste("MESSAGE:", messages),
                   if (length(warnings)) paste("WARNING:", warnings),
                   if (!is.null(failure)) paste("ERROR:", failure)),
                 file.path(out, "logs", paste0(label, ".log")))
    status <- if (!is.null(failure)) "error" else if (length(warnings)) "warning" else "pass"
    record(label, kind, status, start, end, warnings = warnings,
            error = if (is.null(failure)) "" else failure)
    identical(status, "pass")
  }
  run_code <- function(code, environment) {
    expressions <- parse(text = code, keep.source = FALSE)
    # Refuse accidental install/download instructions rather than inferring
    # that a whole chunk should be skipped from a matching comment/string.
    forbidden <- c("install.packages", "download.file", "install_github",
                    "install_git", "install_url", "pkg_install", "curl_download")
    calls <- function(x) {
      if (!is.call(x)) return(character())
      head <- x[[1L]]
      name <- if (is.symbol(head)) as.character(head) else if (is.call(head) &&
        length(head) == 3L && is.symbol(head[[1L]]) &&
        as.character(head[[1L]]) %in% c("::", ":::")) as.character(head[[3L]]) else ""
      c(name, unlist(lapply(as.list(x)[-1L], calls), use.names = FALSE))
    }
    prohibited <- intersect(unique(unlist(lapply(expressions, calls))), forbidden)
    if (length(prohibited)) {
      stop("Install/download instructions must be explicitly marked eval: false: ",
           paste(prohibited, collapse = ", "), call. = FALSE)
    }
    for (expression in expressions) {
      visible <- withVisible(eval(expression, envir = environment))
      if (visible$visible) print(visible$value)
    }
  }

  rd_files <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
  html_links <- c(prcomp = "stats-prcomp.html")
  html_topics <- character()
  for (rd in rd_files) {
    topic <- tools::file_path_sans_ext(basename(rd))
    execute(paste0("aliases-", topic), "help-aliases", function() {
      parsed <- tools::parse_Rd(rd)
      aliases <- vapply(Filter(function(x) identical(attr(x, "Rd_tag"), "\\alias"),
                               parsed), paste, collapse = "", FUN.VALUE = character(1))
      html_links[aliases] <<- paste0(topic, ".html")
      title <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\title"), parsed)
      html_topics[topic] <<- trimws(paste(unlist(title), collapse = ""))
    })
  }
  execute("html-preview-assets", "help-html-assets", function() {
    stylesheet <- file.path(R.home("doc"), "html", "R.css")
    if (!file.copy(stylesheet, file.path(out, "html", "R.css"))) {
      stop("Could not copy the installed R HTML stylesheet.", call. = FALSE)
    }
    html_escape <- function(x) {
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      gsub(">", "&gt;", x, fixed = TRUE)
    }
    entries <- sprintf('<li><a href="%s.html">%s</a>: %s</li>',
                         names(html_topics), html_escape(names(html_topics)),
                         html_escape(html_topics))
    writeLines(c('<!DOCTYPE html>', '<html lang="en"><head>',
                   '<meta charset="UTF-8">',
                   paste0('<title>', html_escape(package), ' help preview</title>'),
                   '<link rel="stylesheet" type="text/css" href="R.css">',
                   '</head><body><main>',
                   paste0('<h1>', html_escape(package), ' ',
                            html_escape(installed_version), ' help preview</h1>'),
                   '<p>Package topics rendered from the source man pages.</p>',
                   '<ul>', entries, '</ul>',
                   '<h2>Explicit external reference</h2>',
                   paste0('<p><a href="stats-prcomp.html">stats::prcomp</a> ',
                            'is rendered from this R installation. Its further ',
                            'help-topic references are plain text; their pages ',
                            'and examples are not collected or executed.</p>'),
                   '</main></body></html>'),
                 file.path(out, "html", "00Index.html"))
    cat("Stylesheet:", stylesheet, "\n")
  })
  execute("html-external-prcomp", "help-html-reference", function() {
    prcomp_rd <- tools::Rd_db("stats")[["prcomp.Rd"]]
    if (is.null(prcomp_rd)) stop("Installed stats::prcomp help was not found.", call. = FALSE)
    tools::Rd2HTML(prcomp_rd, out = file.path(out, "html", "stats-prcomp.html"),
                    package = "stats", no_links = TRUE,
                    stages = c("build", "install", "render"))
    cat("Reference source: installed stats", as.character(utils::packageVersion("stats")),
          "at", find.package("stats"), "\n")
  })
  # Resolve only the declared external stats::prcomp reference. Working on
  # the parsed Rd avoids rewriting generated HTML or crawling other packages.
  local_preview_links <- function(x) {
    if (is.list(x)) {
      if (identical(attr(x, "Rd_tag"), "\\link") &&
          identical(as.character(attr(x, "Rd_option")), "stats") &&
          identical(paste(unlist(x), collapse = ""), "prcomp")) {
        attr(x, "Rd_option") <- NULL
      }
      for (i in seq_along(x)) x[[i]] <- local_preview_links(x[[i]])
    }
    x
  }
  for (rd in rd_files) {
    topic <- tools::file_path_sans_ext(basename(rd))
    label <- paste0("help-", topic)
    code_path <- file.path(out, "code", paste0(label, ".R"))
    environment <- new.env(parent = example_parent)
    execute(label, "help-example", function() {
      tools::Rd2ex(rd, out = code_path)
      if (file.exists(code_path)) {
        run_code(readLines(code_path, warn = FALSE), environment)
      }
    }, plot = TRUE)
    if (!file.exists(code_path) && identical(tail(results, 1L)[[1L]]$status, "pass")) {
      results[[length(results)]]$status <- "no-examples"
    }
    execute(paste0("html-", topic), "help-html", function() {
      tools::Rd2HTML(local_preview_links(tools::parse_Rd(rd)),
                      out = file.path(out, "html", paste0(topic, ".html")),
                      package = package, Links = html_links,
                      stages = c("build", "install", "render"))
    })
  }

  parse_marker <- function(line, number) {
    text <- sub("^\\s*<!--\\s*", "", line, perl = TRUE)
    text <- sub("\\s*-->\\s*$", "", text, perl = TRUE)
    fields <- strsplit(text, ";", fixed = TRUE)[[1L]]
    values <- list()
    for (field in fields) {
      matched <- regmatches(field, regexec("^\\s*([a-z]+)\\s*:\\s*(.*?)\\s*$",
                                             field, perl = TRUE))[[1L]]
      if (length(matched) != 3L || !matched[2L] %in% c("example", "eval", "reason") ||
          matched[2L] %in% names(values)) {
        stop("Invalid or duplicate example field at README line ", number, call. = FALSE)
      }
      values[[matched[2L]]] <- matched[3L]
    }
    if (is.null(values$example) || !grepl("^[a-z][a-z0-9-]*$", values$example)) {
      stop("Use a lowercase hyphenated example ID at README line ", number, call. = FALSE)
    }
    evaluated <- if (is.null(values$eval)) "true" else values$eval
    if (!evaluated %in% c("true", "false")) {
      stop("eval must be true or false at README line ", number, call. = FALSE)
    }
    reason <- if (is.null(values$reason)) "" else values$reason
    if (evaluated == "false" && !nzchar(reason)) {
      if (grepl("^install($|-|ation)", values$example)) {
        reason <- "Installation instructions explicitly reserved for manual execution."
      } else stop("Unevaluated external example needs a reason at README line ",
                  number, call. = FALSE)
    }
    list(id = values$example, evaluate = evaluated == "true", reason = reason,
         marker_line = number)
  }
  parse_readme <- function(lines) {
    blocks <- list()
    pending <- NULL
    fence <- NULL
    ids <- character()
    for (i in seq_along(lines)) {
      line <- lines[i]
      if (!is.null(fence)) {
        closing <- paste0("^\\s{0,3}", strrep(fence$character, fence$width),
                           fence$character, "*\\s*$")
        if (grepl(closing, line, perl = TRUE)) {
          if (fence$is_r) {
            code <- if (i > fence$line + 1L) lines[seq.int(fence$line + 1L, i - 1L)] else character()
            blocks[[length(blocks) + 1L]] <- c(fence$marker,
              list(start = fence$line + 1L, end = i - 1L, code = code))
          }
          fence <- NULL
        }
        next
      }
      if (grepl("<!--\\s*example\\s*:", line, perl = TRUE)) {
        if (!grepl("^\\s*<!--.*-->\\s*$", line, perl = TRUE) || !is.null(pending)) {
          stop("Malformed or orphaned example marker at README line ", i, call. = FALSE)
        }
        pending <- parse_marker(line, i)
        if (pending$id %in% ids) stop("Duplicate README example ID: ", pending$id, call. = FALSE)
        ids <- c(ids, pending$id)
        next
      }
      opening <- regmatches(line, regexec("^\\s{0,3}(`{3,}|~{3,})\\s*(.*?)\\s*$",
                                             line, perl = TRUE))[[1L]]
      if (length(opening)) {
        is_r <- tolower(opening[3L]) == "r"
        if (is_r && is.null(pending)) {
          stop("Unlabeled R chunk at README line ", i,
               "; add <!-- example: unique-id --> immediately before it.", call. = FALSE)
        }
        if (!is_r && !is.null(pending)) {
          stop("Example marker must precede an r fence at README line ", i, call. = FALSE)
        }
        fence <- list(character = substr(opening[2L], 1L, 1L), width = nchar(opening[2L]),
                       line = i, is_r = is_r, marker = pending)
        pending <- NULL
      } else if (!is.null(pending) && nzchar(trimws(line))) {
        stop("Example marker is not followed by a code fence at README line ", i, call. = FALSE)
      }
    }
    if (!is.null(fence)) stop("Unclosed README fence beginning at line ", fence$line, call. = FALSE)
    if (!is.null(pending)) stop("Orphaned example marker at end of README.", call. = FALSE)
    if (!length(blocks)) stop("No labeled README R examples were found.", call. = FALSE)
    blocks
  }
  blocks <- NULL
  execute("README-parser", "README-parser", function() {
    blocks <<- parse_readme(readLines(file.path(root, "README.md"), warn = FALSE))
  })
  readme_environment <- new.env(parent = example_parent)
  blocked <- FALSE
  for (block in blocks) {
    label <- paste0("README-", block$id)
    writeLines(block$code, file.path(out, "code", paste0(label, ".R")))
    if (!block$evaluate) {
      record(label, "README-example", "skip-explicit", block$start, block$end,
             reason = block$reason)
      writeLines(paste("NOT EXECUTED:", block$reason),
                   file.path(out, "logs", paste0(label, ".log")))
    } else if (blocked) {
      record(label, "README-example", "blocked-after-error", block$start, block$end,
             reason = "An earlier sequential README chunk did not pass.")
    } else {
      passed <- execute(label, "README-example", function() {
        run_code(block$code, readme_environment)
      }, block$start, block$end, plot = TRUE)
      if (!passed) blocked <- TRUE
    }
  }
  source_after <- hashes(source_files)
  manifest$source_files_unchanged_during_run <- identical(source_before, source_after)
  if (!manifest$source_files_unchanged_during_run) {
    record("source-consistency", "provenance", "error",
             error = "Source files changed during this run; rerun after editing is complete.")
  }
  manifest$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  manifest$README_chunks <- length(blocks)
  manifest$README_evaluated_chunks <- sum(vapply(blocks, `[[`, logical(1), "evaluate"))
  saveRDS(manifest, file.path(out, "run-manifest.rds"))
  writeLines(capture.output(dput(manifest)), file.path(out, "run-manifest.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
  result <- do.call(rbind, results)
  utils::write.csv(result, file.path(out, "example-results.csv"), row.names = FALSE)
  print(result, row.names = FALSE)
  if (any(result$status %in% c("error", "warning", "blocked-after-error"))) {
    stop("Documentation validation failed; inspect example-results.csv and logs/.", call. = FALSE)
  }
  invisible(result)
}

main()
