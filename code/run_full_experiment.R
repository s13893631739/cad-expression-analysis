script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))

source(file.path(script_dir, "01_check_environment.R"))
cat("\n=== Running 20_full_pipeline.R ===\n")
source(file.path(script_dir, "20_full_pipeline.R"))

experiment_root <- file.path(output_root, "20_full_experiment")
summary_path <- file.path(experiment_root, "summary.txt")
summary <- read.delim(summary_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
historical_intersection <- length(readLines(file.path(project_root, "results", "08_venn", "intersectGenes.txt"), warn = FALSE))
writeLines(
  c(
    "This is a fresh exploratory run from the restored project inputs.",
    paste("Fresh differential-expression count:", summary$value[summary$item == "deg_filtered"]),
    paste("Historical differential-expression count:", "284"),
    paste("Fresh intersection count:", summary$value[summary$item == "intersect_genes"]),
    paste("Historical intersection count:", historical_intersection),
    "Fresh outputs are intentionally not used to overwrite the historical results/ snapshot."
  ),
  file.path(experiment_root, "RECOMPUTATION_NOTE.txt")
)

session_path <- file.path(experiment_root, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_path)
cat("\nFull experimental path finished. Its outputs are in generated/20_full_experiment.\n")
cat("These outputs are intentionally separate from the saved-input evaluation path run_final.R.\n")
