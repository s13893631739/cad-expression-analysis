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

cat("\nFinal analysis pipeline finished.\n")
