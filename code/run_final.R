script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))

scripts <- c(
  "01_check_environment.R",
  "30_generate_figures.R",
  "40_external_validation.R"
)

for (script in scripts) {
  cat("\n=== Running", script, "===\n")
  source(file.path(script_dir, script))
}

output_files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
if (length(output_files) == 0) {
  stop("The final pipeline produced no output files.")
}
writeLines(
  c(
    "This run evaluates the saved intermediate inputs in data/intermediate.",
    "It does not run code/20_full_pipeline.R or claim a fresh end-to-end reconstruction.",
    paste("Output files before metadata:", length(output_files))
  ),
  file.path(output_root, "RUN_METADATA.txt")
)
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))

cat("\nFinal analysis pipeline finished.\n")
