script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(pROC)
})

external_out <- file.path(output_root, "11_external_validation")
if (dir.exists(external_out)) {
  unlink(list.files(external_out, full.names = TRUE, all.files = TRUE, no.. = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(external_out, recursive = TRUE, showWarnings = FALSE)

series_path <- paths$gse113079_series_matrix
platform_path <- paths$gse113079_platform
hub_genes <- read_gene_list(paths$roc_genes)

split_gene_symbols <- function(x) {
  if (is.na(x) || !nzchar(x) || x %in% c("NA", "None")) {
    return(character(0))
  }
  trimws(unlist(strsplit(x, "[,;|/]+")))
}

read_series_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1 || end <= begin + 1) {
    stop("Unable to locate series matrix table section in ", path)
  }

  expr_df <- read.table(
    path,
    skip = begin,
    nrows = end - begin - 2,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  rownames(expr_df) <- expr_df[[1]]
  expr <- data.matrix(expr_df[, -1, drop = FALSE])

  title_idx <- grep("^!Sample_title", lines)
  if (length(title_idx) != 1) {
    stop("Unable to locate sample titles in ", path)
  }
  sample_titles <- gsub("^\"|\"$", "", strsplit(lines[title_idx], "\t", fixed = TRUE)[[1]][-1])
  group <- ifelse(grepl("healthy control", sample_titles, ignore.case = TRUE), "Control", "CAD")

  list(
    expr = expr,
    group = factor(group, levels = c("Control", "CAD")),
    sample_titles = sample_titles,
    raw_lines = lines
  )
}

read_platform <- function(path) {
  read.delim(
    path,
    skip = 41,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

collapse_gene_expression <- function(expr, platform, gene, agg = c("mean", "median")) {
  agg <- match.arg(agg)
  hits <- vapply(platform$GeneSymbol, function(x) gene %in% split_gene_symbols(x), logical(1))
  probe_ids <- unique(na.omit(platform$ID[hits]))
  probe_ids <- intersect(probe_ids, rownames(expr))
  if (length(probe_ids) == 0) {
    stop("No probe mapped to gene: ", gene)
  }

  vals <- expr[probe_ids, , drop = FALSE]
  if (length(probe_ids) == 1) {
    vec <- as.numeric(vals[1, ])
  } else if (agg == "mean") {
    vec <- colMeans(vals, na.rm = TRUE)
  } else {
    vec <- apply(vals, 2, median, na.rm = TRUE)
  }
  list(values = vec, probes = probe_ids)
}

save_base_plot <- function(path_no_ext, width = 7, height = 7, draw_fun, res = 300) {
  pdf(paste0(path_no_ext, ".pdf"), width = width, height = height, useDingbats = FALSE)
  draw_fun()
  dev.off()

  png(paste0(path_no_ext, ".png"), width = width * res, height = height * res, res = res)
  draw_fun()
  dev.off()
}

roc_xy <- function(roc_obj) {
  xy <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  xy <- xy[order(xy$fpr, xy$tpr), , drop = FALSE]
  xy <- unique(rbind(data.frame(fpr = 0, tpr = 0), xy, data.frame(fpr = 1, tpr = 1)))
  xy
}

draw_step_roc <- function(roc_obj, col, lwd = 2.8) {
  xy <- roc_xy(roc_obj)
  lines(xy$fpr, xy$tpr, type = "s", col = col, lwd = lwd)
}

external <- read_series_matrix(series_path)
platform <- read_platform(platform_path)
expr <- external$expr
group <- external$group

gene_expr <- matrix(NA_real_, nrow = ncol(expr), ncol = length(hub_genes))
colnames(gene_expr) <- hub_genes
rownames(gene_expr) <- colnames(expr)

probe_map <- data.frame(
  Gene = character(),
  Probes = character(),
  nProbes = integer(),
  Aggregation = character(),
  stringsAsFactors = FALSE
)

for (gene in hub_genes) {
  mapped <- collapse_gene_expression(expr, platform, gene, agg = "mean")
  gene_expr[, gene] <- mapped$values
  probe_map <- rbind(
    probe_map,
    data.frame(
      Gene = gene,
      Probes = paste(mapped$probes, collapse = ";"),
      nProbes = length(mapped$probes),
      Aggregation = "mean",
      stringsAsFactors = FALSE
    )
  )
}

roc_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)

train_expr <- read.table(paths$merged_after_combat, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
train_meta <- read.table(paths$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
train_group <- factor(train_meta$group, levels = c("Control", "Treat"))

train_panel <- data.frame(
  group = train_group,
  t(train_expr[hub_genes, train_meta$sample, drop = FALSE]),
  check.names = FALSE
)

roc_objects <- list()
auc_table <- data.frame(
  Gene = character(),
  AUC = numeric(),
  nProbes = integer(),
  Probes = character(),
  stringsAsFactors = FALSE
)

for (gene in hub_genes) {
  train_direction <- if (
    mean(train_panel[train_panel$group == "Treat", gene], na.rm = TRUE) >=
      mean(train_panel[train_panel$group == "Control", gene], na.rm = TRUE)
  ) "<" else ">"
  roc_obj <- roc(
    response = group,
    predictor = gene_expr[, gene],
    levels = c("Control", "CAD"),
    direction = train_direction,
    quiet = TRUE
  )
  roc_objects[[gene]] <- roc_obj
  auc_table <- rbind(
    auc_table,
    data.frame(
      Gene = gene,
      AUC = as.numeric(auc(roc_obj)),
      Direction = train_direction,
      nProbes = probe_map$nProbes[probe_map$Gene == gene],
      Probes = probe_map$Probes[probe_map$Gene == gene],
      stringsAsFactors = FALSE
    )
  )
}

plot_external_roc <- function() {
  par(mar = c(5.2, 5.2, 4.2, 2.0), xaxs = "i", yaxs = "i")
  plot(
    NA,
    main = "External Validation ROC Curves of 10 Hub Genes in GSE113079",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    xlim = c(0, 1),
    ylim = c(0, 1),
    las = 1
  )
  for (i in seq_along(hub_genes)) {
    draw_step_roc(roc_objects[[hub_genes[i]]], col = roc_colors[i], lwd = 2.8)
  }
  abline(a = 0, b = 1, lty = 2, col = "#377eb8", lwd = 1.6)
  legend_labels <- sprintf("%s (AUC=%.3f)", auc_table$Gene, auc_table$AUC)
  legend(
    "bottomright",
    legend = legend_labels,
    col = roc_colors,
    lwd = 2.8,
    cex = 0.95,
    bty = "n"
  )
}

save_base_plot(
  file.path(external_out, "GSE113079_external_ROC_10_hub_genes"),
  width = 8,
  height = 8,
  draw_fun = plot_external_roc,
  res = 300
)

panel_fit <- glm(group ~ ., data = train_panel, family = binomial())
panel_coefs <- coef(panel_fit)

# Match the raw feature scale used when fitting panel_fit.
# Cross-platform comparability still needs evaluation; no historical AUC is reused.
panel_logit_score <- as.numeric(predict(
  panel_fit, newdata = as.data.frame(gene_expr[, hub_genes, drop = FALSE]),
  type = "link"
))
panel_train_score <- as.numeric(predict(panel_fit, newdata = train_panel, type = "link"))
panel_train_direction <- if (
  mean(panel_train_score[train_group == "Treat"], na.rm = TRUE) >=
    mean(panel_train_score[train_group == "Control"], na.rm = TRUE)
) "<" else ">"
panel_logit_roc <- roc(group, panel_logit_score, levels = c("Control", "CAD"), direction = panel_train_direction, quiet = TRUE)
panel_logit_auc <- as.numeric(auc(panel_logit_roc))

train_direction <- vapply(hub_genes, function(g) {
  if (
    mean(train_panel[train_panel$group == "Treat", g], na.rm = TRUE) >=
      mean(train_panel[train_panel$group == "Control", g], na.rm = TRUE)
  ) 1 else -1
}, numeric(1))

internal_auc <- read.table(
  file.path(output_root, "10_roc", "roc_auc_table.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
internal_weights <- setNames(internal_auc$AUC - 0.5, internal_auc$Gene)
internal_weights <- internal_weights[hub_genes]

train_x <- as.matrix(train_panel[, hub_genes, drop = FALSE])
train_center <- colMeans(train_x, na.rm = TRUE)
train_scale <- apply(train_x, 2, sd, na.rm = TRUE)
train_scale[!is.finite(train_scale) | train_scale == 0] <- 1
ext_scaled_train <- scale(gene_expr[, hub_genes, drop = FALSE], center = train_center, scale = train_scale)
train_scaled <- scale(train_x, center = train_center, scale = train_scale)
weighted_train_score <- as.numeric(train_scaled %*% (train_direction * internal_weights))
weighted_score <- as.numeric(ext_scaled_train %*% (train_direction * internal_weights))
weighted_score <- weighted_score / sum(abs(internal_weights))
weighted_train_score <- weighted_train_score / sum(abs(internal_weights))
weighted_train_direction <- if (
  mean(weighted_train_score[train_group == "Treat"], na.rm = TRUE) >=
    mean(weighted_train_score[train_group == "Control"], na.rm = TRUE)
) "<" else ">"
weighted_roc <- roc(group, weighted_score, levels = c("Control", "CAD"), direction = weighted_train_direction, quiet = TRUE)
weighted_auc <- as.numeric(auc(weighted_roc))

plot_panel_roc <- function() {
  par(mar = c(5.2, 5.2, 4.2, 2.0), xaxs = "i", yaxs = "i")
  plot(
    NA,
    main = "External Validation ROC of 10-Gene Signature Score",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    xlim = c(0, 1),
    ylim = c(0, 1),
    las = 1
  )
  draw_step_roc(weighted_roc, col = "#1f77b4", lwd = 2.8)
  draw_step_roc(panel_logit_roc, col = "#d62728", lwd = 2.4)
  abline(a = 0, b = 1, lty = 2, col = "#377eb8", lwd = 1.6)
  legend(
    "bottomright",
    legend = c(
      sprintf("Weighted score (AUC=%.3f)", weighted_auc),
      sprintf("Transferred GLM (AUC=%.3f)", panel_logit_auc)
    ),
    col = c("#1f77b4", "#d62728"),
    lwd = 2.8,
    cex = 0.95,
    bty = "n"
  )
}

save_base_plot(
  file.path(external_out, "GSE113079_external_ROC_10_gene_signature"),
  width = 7.5,
  height = 7.5,
  draw_fun = plot_panel_roc,
  res = 300
)

summary_tbl <- data.frame(
  item = c(
    "dataset", "samples_total", "samples_cad", "samples_control",
    "single_gene_reference_fig", "weighted_signature_auc", "transfer_glm_auc",
    "weighted_score_scaling", "weighted_score_direction", "transfer_glm_direction"
  ),
  value = c(
    "GSE113079", ncol(expr), sum(group == "CAD"), sum(group == "Control"),
    "GSE113079_external_ROC_10_hub_genes", sprintf("%.3f", weighted_auc), sprintf("%.3f", panel_logit_auc),
    "training_cohort_center_scale", weighted_train_direction, panel_train_direction
  ),
  stringsAsFactors = FALSE
)

write.table(
  auc_table,
  file.path(external_out, "roc_auc_table.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  probe_map,
  file.path(external_out, "probe_mapping.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_tbl,
  file.path(external_out, "external_validation_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- data.frame(
  file = c(
    "GSE113079_external_ROC_10_hub_genes.pdf",
    "GSE113079_external_ROC_10_hub_genes.png",
    "GSE113079_external_ROC_10_gene_signature.pdf",
    "GSE113079_external_ROC_10_gene_signature.png",
    "roc_auc_table.txt",
    "probe_mapping.txt",
    "external_validation_summary.txt"
  ),
  description = c(
    "Ten individual hub-gene ROC curves on GSE113079",
    "Ten individual hub-gene ROC curves on GSE113079",
    "10-gene weighted signature ROC on GSE113079",
    "10-gene weighted signature ROC on GSE113079",
    "Single-gene AUC summary",
    "Probe-to-gene aggregation map",
    "Run summary"
  ),
  stringsAsFactors = FALSE
)

write.table(
  manifest,
  file.path(external_out, "validation_manifest.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("External validation finished. Outputs:", external_out, "\n")
