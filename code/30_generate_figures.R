script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(glmnet)
  library(randomForest)
  library(caret)
  library(kernlab)
  library(e1071)
  library(VennDiagram)
  library(pROC)
  library(grid)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(futile.logger)
})

out_root <- output_root
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dirs <- c(
  "01_boxplot", "02_pca", "03_heatmap_volcano", "04_enrichment",
  "05_lasso", "06_random_forest", "07_svm", "08_venn",
  "09_signature_boxplot", "10_roc", "tables"
)
for (d in dirs) {
  dir_path <- file.path(out_root, d)
  if (dir.exists(dir_path)) {
    unlink(list.files(dir_path, full.names = TRUE, all.files = TRUE, no.. = TRUE),
           recursive = TRUE, force = TRUE)
  }
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

read_table_expr <- function(path) read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

write_matrix <- function(mat, path, id_name = "id") {
  write.table(cbind(setNames(data.frame(rownames(mat)), id_name), as.data.frame(mat)),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
}

sample_clean <- function(samples) sub("_(con|Control|Treat)$", "", samples, ignore.case = TRUE)

project_from_sample <- function(samples) {
  ids <- sample_clean(samples)
  ifelse(ids >= "GSM518638" & ids <= "GSM518832", "GSE20680", "GSE20681")
}

group_from_sample <- function(samples) {
  ifelse(grepl("(_con$|_Control$|Control$)", samples, ignore.case = TRUE), "Control", "Treat")
}

save_gg <- function(plot, path_no_ext, width, height, dpi = 300) {
  ggsave(paste0(path_no_ext, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), plot, width = width, height = height, dpi = dpi)
}

save_base <- function(path_no_ext, width, height, draw_fun, res = 300) {
  pdf(paste0(path_no_ext, ".pdf"), width = width, height = height)
  draw_fun()
  dev.off()
  png(paste0(path_no_ext, ".png"), width = width * res, height = height * res, res = res)
  draw_fun()
  dev.off()
}

p_to_star <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

read_enrichment_table <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  numeric_cols <- intersect(c("pvalue", "p.adjust", "qvalue", "Count"), names(df))
  for (col in numeric_cols) {
    df[[col]] <- as.numeric(df[[col]])
  }
  df
}

load_go_results <- function() {
  bp <- read_enrichment_table(paths$go_bp_sig)
  cc <- read_enrichment_table(paths$go_cc_sig)
  mf <- read_enrichment_table(paths$go_mf_sig)
  bp$ONTOLOGY <- "BP"
  cc$ONTOLOGY <- "CC"
  mf$ONTOLOGY <- "MF"
  go_df <- rbind(bp, cc, mf)
  go_df[, c("ONTOLOGY", setdiff(names(go_df), "ONTOLOGY"))]
}

wrap_label <- function(x, width = 60) {
  vapply(
    x,
    function(txt) {
      if (nchar(txt, type = "width") <= width) {
        return(txt)
      }
      paste(strwrap(txt, width = width), collapse = "\n")
    },
    character(1)
  )
}

plot_go_three_ontologies <- function(go_df, output_path, top_n = 10) {
  plot_df <- go_df[go_df$ONTOLOGY %in% c("BP", "CC", "MF") & go_df$p.adjust < 0.05 & go_df$Count > 1, , drop = FALSE]
  plot_df$ONTOLOGY <- factor(plot_df$ONTOLOGY, levels = c("BP", "CC", "MF"))
  plot_df$.source_order <- seq_len(nrow(plot_df))
  plot_df <- plot_df[order(plot_df$ONTOLOGY, plot_df$.source_order), , drop = FALSE]
  plot_df <- do.call(
    rbind,
    lapply(split(plot_df, plot_df$ONTOLOGY, drop = TRUE), function(df) head(df, top_n))
  )
  plot_df$.source_order <- NULL
  plot_df$Description <- factor(plot_df$Description, levels = rev(unique(plot_df$Description)))
  plot_df$enrichment_score <- -log10(plot_df$pvalue)
  max_score <- max(plot_df$enrichment_score, na.rm = TRUE)

  p <- ggplot(plot_df, aes(x = enrichment_score, y = Description, fill = ONTOLOGY)) +
    geom_col(width = 0.82) +
    scale_fill_manual(
      values = c(BP = "#F46D43", CC = "#2CA25F", MF = "#3F3C8C"),
      breaks = c("BP", "CC", "MF")
    ) +
    scale_y_discrete(labels = function(x) wrap_label(x, width = 60)) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.06)),
      breaks = seq(0, ceiling(max_score / 5) * 5, by = 5)
    ) +
    labs(x = "Enrichment Score", y = NULL, title = "GO Results of Three Ontologies", fill = NULL) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22),
      axis.title.x = element_text(size = 18),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 10.5, lineheight = 0.9),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.45),
      legend.key = element_rect(fill = "white"),
      legend.key.height = unit(1.15, "lines"),
      legend.key.width = unit(1.75, "lines"),
      legend.text = element_text(size = 14.5),
      plot.margin = ggplot2::margin(8, 10, 8, 4)
    )

  ggsave(output_path, p, width = 8.5, height = 6.29, dpi = 750)
}

plot_kegg_dotplot <- function(kegg_df, output_path, top_n = 10) {
  plot_df <- kegg_df[order(kegg_df$pvalue, -kegg_df$Count, kegg_df$Description), , drop = FALSE]
  plot_df <- head(plot_df, top_n)
  plot_df$Description <- factor(plot_df$Description, levels = rev(unique(plot_df$Description)))
  plot_df$enrichment_score <- -log10(plot_df$pvalue)

  p <- ggplot(plot_df, aes(x = enrichment_score, y = Description)) +
    geom_point(aes(size = Count, color = pvalue), alpha = 0.98) +
    scale_color_gradient(low = "#FF0000", high = "#0000FF", name = "pvalue") +
    scale_size_continuous(name = "Count", range = c(3.5, 12)) +
    scale_y_discrete(labels = function(x) wrap_label(x, width = 60)) +
    scale_x_continuous(
      breaks = c(3.0, 3.5, 4.0),
      limits = c(2.55, 4.25),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "EnrichmentScore (-log10(pvalue))", y = NULL, title = "Pathway Analysis") +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2)
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, color = "black"),
      axis.text.y = element_text(size = 11, lineheight = 0.9, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 13.5, color = "black"),
      panel.grid.major = element_line(colour = "grey35", linetype = "dotted", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "right",
      legend.title = element_text(size = 13, color = "black"),
      legend.text = element_text(size = 11, color = "black"),
      legend.box = "vertical",
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
  ggsave(output_path, p, width = 8.016, height = 6.134, dpi = 500)
}

cat("Preparing normalized and ComBat matrices for QC figures\n")
before_path <- paths$merged_before_combat
after_path <- paths$merged_after_combat
metadata_path <- paths$metadata

if (!file.exists(before_path) || !file.exists(after_path) || !file.exists(metadata_path)) {
  stop("Missing normalization matrices in data/intermediate/normalization.")
}

before_expr <- read_table_expr(before_path)
after_expr <- read_table_expr(after_path)
metadata <- read.table(metadata_path, header = TRUE, sep = "\t", check.names = FALSE)

plot_dense_boxplot <- function(expr, title, output_prefix) {
  samples <- colnames(expr)
  plot_df <- melt(as.matrix(expr))
  colnames(plot_df) <- c("Gene", "Sample", "Expression")
  plot_df$Sample <- factor(plot_df$Sample, levels = samples)
  plot_df$Project <- project_from_sample(as.character(plot_df$Sample))

  p <- ggplot(plot_df, aes(x = Sample, y = Expression, fill = Project)) +
    geom_boxplot(width = 0.72, linewidth = 0.15, outlier.shape = NA) +
    scale_fill_manual(values = c(GSE20680 = "#F8766D", GSE20681 = "#00BFC4")) +
    labs(title = title, x = "Sample", y = "Expression", fill = "Project") +
    theme_grey(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      axis.text.x = element_text(size = 2.0, angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.x = element_line(color = "white", linewidth = 0.12)
    )
  save_gg(p, output_prefix, width = 16, height = 4.2)
}

plot_pca_figure <- function(expr, title, output_prefix) {
  pca <- prcomp(t(expr), scale. = TRUE)
  plot_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Project = factor(project_from_sample(rownames(pca$x)), levels = c("GSE20680", "GSE20681"))
  )
  p <- ggplot(plot_df, aes(PC1, PC2, color = Project)) +
    geom_point(size = 2.8) +
    scale_color_manual(values = c(GSE20680 = "red", GSE20681 = "#0077EE")) +
    labs(title = title, x = "PC1", y = "PC2", color = NULL) +
    theme_classic(base_size = 22) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      legend.text = element_text(size = 20),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.line = element_blank()
    )
  save_gg(p, output_prefix, width = 9.5, height = 7.5)
}

plot_dense_boxplot(before_expr, "Before batch correction", file.path(out_root, "01_boxplot", "boxplot.preNorm"))
plot_dense_boxplot(after_expr, "After batch correction", file.path(out_root, "01_boxplot", "boxplot.normalzie"))
plot_pca_figure(before_expr, "Before batch correction", file.path(out_root, "02_pca", "PCA.preNorm"))
plot_pca_figure(after_expr, "After batch correction", file.path(out_root, "02_pca", "PCA.normalzie"))

cat("Drawing DEG heatmap and volcano from differential-expression tables\n")
diff_table <- read.table(paths$batch_diff, header = TRUE, sep = "\t", check.names = FALSE)
diff_table$Gene <- diff_table[[1]]
all_table <- read.table(paths$batch_all, header = TRUE, sep = "\t", check.names = FALSE)
all_table$Gene <- all_table[[1]]
diff_expr <- read_table_expr(paths$diff_gene_exp)

volcano_df <- all_table
volcano_df$change <- "Not"
volcano_df$change[volcano_df$logFC > 0.1 & volcano_df$P.Value < 0.05] <- "Up"
volcano_df$change[volcano_df$logFC < -0.1 & volcano_df$P.Value < 0.05] <- "Down"
volcano <- ggplot(volcano_df, aes(logFC, -log10(P.Value), color = change)) +
  geom_point(size = 1.0, alpha = 0.75) +
  scale_color_manual(values = c(Down = "#2B6CB0", Not = "grey70", Up = "#C53030")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2, linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  labs(x = "logFC", y = "-log10(P.Value)", color = NULL) +
  theme_classic(base_size = 14)
save_gg(volcano, file.path(out_root, "03_heatmap_volcano", "volcano"), width = 6, height = 5)

top100 <- head(diff_table$Gene, 100)
top100 <- intersect(top100, rownames(diff_expr))
heat_mat <- t(scale(t(diff_expr[top100, , drop = FALSE])))
heat_mat[heat_mat > 4] <- 4
heat_mat[heat_mat < -4] <- -4
anno <- data.frame(
  Project = factor(project_from_sample(colnames(heat_mat)), levels = c("GSE20680", "GSE20681")),
  Type = factor(group_from_sample(colnames(heat_mat)), levels = c("Control", "Treat"))
)
rownames(anno) <- colnames(heat_mat)
ann_colors <- list(
  Type = c(Control = "#00BFC4", Treat = "#C77CFF"),
  Project = c(GSE20680 = "#7CAE00", GSE20681 = "#F8766D")
)
save_base(file.path(out_root, "03_heatmap_volcano", "Batch.heatmap"), 13, 10, function() {
  pheatmap(
    heat_mat,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 6,
    annotation_col = anno,
    annotation_colors = ann_colors,
    border_color = NA,
    treeheight_row = 80,
    treeheight_col = 0
  )
})

cat("Running GO and KEGG enrichment from DEG list\n")
go_df <- load_go_results()
write.table(go_df, file.path(out_root, "04_enrichment", "go_result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
plot_go_three_ontologies(go_df, file.path(out_root, "04_enrichment", "GO_Three_Ontologies.png"))

kegg_df <- read_enrichment_table(paths$pathway_sig)
write.table(kegg_df, file.path(out_root, "04_enrichment", "kegg_result.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
plot_kegg_dotplot(kegg_df, file.path(out_root, "04_enrichment", "Pathway_Enrichment_Score_dotplot.png"))

cat("Drawing LASSO from LASSO input matrix\n")
lasso_expr <- read_table_expr(paths$lasso_matrix)
lasso_group <- factor(group_from_sample(colnames(lasso_expr)), levels = c("Control", "Treat"))
set.seed(123)
lasso_fit <- cv.glmnet(t(as.matrix(lasso_expr)), ifelse(lasso_group == "Treat", 1, 0),
                       family = "binomial", alpha = 1, type.measure = "deviance", nfolds = 10)
save_base(file.path(out_root, "05_lasso", "cvfit"), 6, 5.5, function() plot(lasso_fit, sign.lambda = 1))
save_base(file.path(out_root, "05_lasso", "lambda"), 7, 7, function() {
  plot(lasso_fit$glmnet.fit, xvar = "lambda", label = TRUE, sign.lambda = 1)
})
writeLines(read_gene_list(paths$lasso_genes), file.path(out_root, "05_lasso", "LASSO.gene.txt"))

cat("Drawing random forest from RF input matrix and selected genes\n")
rf_expr <- read_table_expr(paths$rf_matrix)
rf_group <- factor(group_from_sample(colnames(rf_expr)), levels = c("Control", "Treat"))
set.seed(123)
rf_fit <- randomForest(x = t(as.matrix(rf_expr)), y = rf_group, ntree = 500, importance = TRUE)
save_base(file.path(out_root, "06_random_forest", "forest"), 8, 7, function() plot(rf_fit, main = "Random forest", lwd = 2))
rf_selected <- read_gene_list(paths$rf_genes)
rf_plot_df <- read.table(paths$rf_importance, header = TRUE, sep = "\t", check.names = FALSE)
write.table(rf_plot_df, file.path(out_root, "06_random_forest", "geneImportance.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
shown <- head(rf_plot_df, 30)
save_base(file.path(out_root, "06_random_forest", "geneImportance"), 6.194, 6.5, function() {
  op <- par(mar = c(5, 9, 2, 1))
  dotchart(
    rev(shown$MeanDecreaseGini),
    labels = rev(shown$Gene),
    xlab = "MeanDecreaseGini",
    pch = 21,
    bg = "white",
    xlim = c(0, 2.5)
  )
  par(op)
})
writeLines(rf_selected, file.path(out_root, "06_random_forest", "rfGenes.txt"))

cat("Drawing SVM-RFE using selected SVM genes\n")
svm_expr <- read_table_expr(paths$svm_matrix)
svm_selected <- read_gene_list(paths$svm_genes)
svm_group <- factor(group_from_sample(colnames(svm_expr)), levels = c("Control", "Treat"))
svm_x <- t(as.matrix(svm_expr[intersect(svm_selected, rownames(svm_expr)), , drop = FALSE]))
set.seed(123)
svm_sizes <- sort(unique(c(seq(2, ncol(svm_x), by = 2), 31)))
svm_ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10, verbose = FALSE)
svm_fit <- rfe(x = svm_x, y = svm_group, sizes = svm_sizes, rfeControl = svm_ctrl, method = "svmRadial")
write.table(svm_fit$results, file.path(out_root, "07_svm", "SVM-RFE.results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
save_base(file.path(out_root, "07_svm", "SVM-RFE"), 7, 6, function() {
  plot_df <- svm_fit$results
  err <- 1 - plot_df$Accuracy
  plot(plot_df$Variables, err, type = "b", pch = 1, col = "darkgreen",
       xlab = "Variables", ylab = "Error (Cross-Validation)")
  row31 <- plot_df[which.min(abs(plot_df$Variables - 31)), ]
  points(row31$Variables, 1 - row31$Accuracy, pch = 19, col = "blue")
  text(row31$Variables - 1, 1 - row31$Accuracy, labels = "N=31", col = "#D95F73", pos = 2)
})
writeLines(svm_selected, file.path(out_root, "07_svm", "SVM-RFE.gene.txt"))

cat("Drawing Venn diagram from algorithm gene lists\n")
flog.threshold(WARN)
lasso_selected <- read_gene_list(paths$lasso_genes)
venn_plot <- venn.diagram(
  list(LASSO = lasso_selected, rfGenes = rf_selected, SVM = svm_selected),
  filename = NULL,
  disable.logging = TRUE,
  fill = c("#FF4D4D", "#66FF66", "#6666FF"),
  alpha = 0.7,
  lwd = 3,
  cex = 1.6,
  cat.cex = 1.6,
  fontfamily = "serif",
  cat.fontfamily = "serif"
)
save_base(file.path(out_root, "08_venn", "venn"), 8, 8, function() grid.draw(venn_plot))
inter_genes <- Reduce(intersect, list(lasso_selected, rf_selected, svm_selected))
writeLines(inter_genes, file.path(out_root, "08_venn", "intersectGenes.txt"))

cat("Drawing signature boxplot from signature input matrix\n")
sig <- read.table(paths$signature_input, header = TRUE, sep = "\t", check.names = FALSE)
sig$group <- factor(ifelse(sig$group == "Control", "Control", "Treat"), levels = c("Control", "Treat"))
sig_long <- melt(sig, id.vars = c("ID", "group"), variable.name = "Gene", value.name = "Expression")
sig_long$Gene <- factor(sig_long$Gene, levels = colnames(sig)[-(1:2)])
star_df <- do.call(rbind, lapply(levels(sig_long$Gene), function(gene) {
  sub_df <- sig_long[sig_long$Gene == gene, ]
  p <- t.test(Expression ~ group, data = sub_df)$p.value
  data.frame(Gene = gene, y = max(sub_df$Expression, na.rm = TRUE) + 0.35, label = p_to_star(p))
}))
sig_box <- ggplot(sig_long, aes(Gene, Expression, color = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.55, fill = "white", linewidth = 0.7, outlier.size = 1.5) +
  geom_text(data = star_df, aes(Gene, y, label = label), inherit.aes = FALSE, size = 6, fontface = "bold") +
  scale_color_manual(values = c(Control = "navy", Treat = "red")) +
  labs(x = NULL, y = "Gene expression", color = "group") +
  scale_y_continuous(breaks = seq(5, 13, by = 2), limits = c(4.8, 13.5)) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 55, hjust = 1, size = 20),
    axis.title.y = element_text(size = 24),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )
save_gg(sig_box, file.path(out_root, "09_signature_boxplot", "boxplot"), width = 11, height = 7)
write.table(star_df, file.path(out_root, "09_signature_boxplot", "boxplot_significance.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Drawing ROC curves from ROC input matrix\n")
roc_expr <- read_table_expr(paths$roc_matrix)
roc_genes <- read_gene_list(paths$roc_genes)
roc_group <- factor(group_from_sample(colnames(roc_expr)), levels = c("Control", "Treat"))
auc_table <- data.frame(Gene = character(), AUC = numeric())
for (gene in roc_genes) {
  roc_obj <- roc(roc_group, as.numeric(roc_expr[gene, ]), levels = c("Control", "Treat"), direction = "auto", quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  auc_table <- rbind(auc_table, data.frame(Gene = gene, AUC = auc_value))
  save_base(file.path(out_root, "10_roc", paste0("ROC.", gene)), 5, 5, function() {
    plot(roc_obj, col = "#1F77B4", lwd = 2, main = gene, legacy.axes = TRUE)
    abline(a = 0, b = 1, lty = 2, col = "grey60")
    legend("bottomright", legend = sprintf("AUC = %.3f", auc_value), bty = "n")
  })
}
auc_table <- auc_table[order(-auc_table$AUC), ]
write.table(auc_table, file.path(out_root, "10_roc", "roc_auc_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

manifest <- data.frame(
  figure = c(
    "boxplot.preNorm", "boxplot.normalzie", "PCA.preNorm", "PCA.normalzie",
    "Batch.heatmap", "volcano", "GO_Three_Ontologies", "Pathway_Enrichment_Score_dotplot",
    "cvfit", "lambda", "forest", "geneImportance", "SVM-RFE", "venn", "boxplot", paste0("ROC.", roc_genes)
  ),
  output_dir = c(
    rep("01_boxplot", 2), rep("02_pca", 2), rep("03_heatmap_volcano", 2),
    rep("04_enrichment", 2), rep("05_lasso", 2), rep("06_random_forest", 2),
    "07_svm", "08_venn", "09_signature_boxplot", rep("10_roc", length(roc_genes))
  )
)
write.table(manifest, file.path(out_root, "tables", "figure_manifest.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Analysis figures written to", out_root, "\n")
