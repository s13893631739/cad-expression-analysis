# Verify that every large input listed in data/input_manifest.tsv is hydrated and intact.
# This script uses base R only and does not train a model.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) == 0) stop("Run this script with Rscript scripts/check_inputs.R from the repository root.")
script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
manifest_path <- file.path(project_root, "data", "input_manifest.tsv")
manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("relative_path", "bytes", "sha256")
if (!all(required_columns %in% names(manifest))) {
  stop("Manifest is missing required columns: ", paste(setdiff(required_columns, names(manifest)), collapse = ", "))
}
if (anyDuplicated(manifest$relative_path)) stop("Manifest contains duplicate relative paths.")

is_lfs_pointer <- function(path) {
  if (!file.exists(path)) return(FALSE)
  first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character())
  length(first_line) == 1 && identical(first_line, "version https://git-lfs.github.com/spec/v1")
}

problems <- character()
for (i in seq_len(nrow(manifest))) {
  relative_path <- manifest$relative_path[i]
  path <- file.path(project_root, relative_path)
  if (!file.exists(path)) {
    problems <- c(problems, paste(relative_path, "is missing"))
    next
  }
  if (is_lfs_pointer(path)) {
    problems <- c(problems, paste(relative_path, "is still a Git LFS pointer; run `git lfs pull`"))
    next
  }
  actual_bytes <- unname(file.info(path)$size)
  expected_bytes <- as.numeric(manifest$bytes[i])
  if (!isTRUE(actual_bytes == expected_bytes)) {
    problems <- c(problems, sprintf("%s has %s bytes; expected %s", relative_path, actual_bytes, expected_bytes))
    next
  }
  actual_sha <- unname(tools::sha256sum(path))
  expected_sha <- manifest$sha256[i]
  if (!identical(tolower(actual_sha), tolower(expected_sha))) {
    problems <- c(problems, paste(relative_path, "has SHA-256", actual_sha, "; expected", expected_sha))
  }
}

cat("Manifest:", manifest_path, "\n")
cat("Inputs checked:", nrow(manifest), "\n")
if (length(problems) > 0) {
  cat("Input verification failed:\n")
  cat(paste0(" - ", problems, collapse = "\n"), "\n", sep = "")
  stop("Fix the input files or run `git lfs pull`, then retry.")
}
cat("All manifest inputs are present, hydrated, and match their byte sizes and SHA-256 checksums.\n")
