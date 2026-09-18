script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(limma)
  library(sva)
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

out_root <- make_out_dir("20_full_experiment")
dir.create(file.path(out_root, "01_normalize"), showWarnings = FALSE)
dir.create(file.path(out_root, "02_qc"), showWarnings = FALSE)
dir.create(file.path(out_root, "03_diff"), showWarnings = FALSE)
dir.create(file.path(out_root, "04_enrichment"), showWarnings = FALSE)
dir.create(file.path(out_root, "05_lasso"), showWarnings = FALSE)
dir.create(file.path(out_root, "06_random_forest"), showWarnings = FALSE)
dir.create(file.path(out_root, "07_svm_rfe"), showWarnings = FALSE)
dir.create(file.path(out_root, "08_intersection"), showWarnings = FALSE)
dir.create(file.path(out_root, "09_signature"), showWarnings = FALSE)
dir.create(file.path(out_root, "10_roc"), showWarnings = FALSE)

write_matrix <- function(mat, path, id_name = "id") {
  write.table(cbind(setNames(data.frame(rownames(mat)), id_name), as.data.frame(mat)),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_samples <- function(path) read_gene_list(path)

make_metadata <- function(samples, dataset, control_samples, treat_samples) {
  data.frame(
    sample = samples,
    clean_sample = clean_sample_id(samples),
    dataset = dataset,
    group = ifelse(clean_sample_id(samples) %in% control_samples, "Control",
                   ifelse(clean_sample_id(samples) %in% treat_samples, "Treat", NA)),
    stringsAsFactors = FALSE
  )
}

project_from_clean_id <- function(ids) {
  ifelse(ids >= "GSM518638" & ids <= "GSM518832", "GSE20680", "GSE20681")
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

plot_dense_boxplot <- function(expr, meta, title, out_prefix) {
  plot_df <- melt(as.matrix(expr))
  colnames(plot_df) <- c("Gene", "Sample", "Expression")
  plot_df$Sample <- factor(plot_df$Sample, levels = colnames(expr))
  plot_df$Project <- meta$dataset[match(as.character(plot_df$Sample), meta$sample)]
  p <- ggplot(plot_df, aes(Sample, Expression, fill = Project)) +
    geom_boxplot(width = 0.72, linewidth = 0.16, outlier.shape = NA) +
    scale_fill_manual(values = c(GSE20680 = "#F8766D", GSE20681 = "#00BFC4")) +
    labs(title = title, x = "Sample", y = "Expression", fill = "Project") +
    theme_grey(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      axis.text.x = element_text(size = 2, angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = element_line(color = "white", linewidth = 0.12),
      legend.position = "right"
    )
  ggsave(paste0(out_prefix, ".pdf"), p, width = 16, height = 4.2)
  ggsave(paste0(out_prefix, ".png"), p, width = 16, height = 4.2, dpi = 300)
}

plot_pca <- function(expr, meta, title, out_prefix) {
  pca <- prcomp(t(expr), scale. = TRUE)
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Project = factor(meta$dataset[match(rownames(pca$x), meta$sample)], levels = c("GSE20680", "GSE20681"))
  )
  p <- ggplot(pca_df, aes(PC1, PC2, color = Project)) +
    geom_point(size = 2.7) +
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
  ggsave(paste0(out_prefix, ".pdf"), p, width = 9.5, height = 7.5)
  ggsave(paste0(out_prefix, ".png"), p, width = 9.5, height = 7.5, dpi = 300)
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
      axis.text.y = element_text(size = 11.5),
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

cat("Step 1: read source matrices and sample groups\n")
raw1 <- read_expression(paths$gse20680_raw)
raw2 <- read_expression(paths$gse20681_raw)
gse20680_con <- read_samples(paths$gse20680_con)
gse20680_treat <- read_samples(paths$gse20680_cad)
gse20681_con <- read_samples(paths$gse20681_con)
gse20681_treat <- read_samples(paths$gse20681_cad)

raw1 <- raw1[, c(gse20680_con, gse20680_treat), drop = FALSE]
raw2 <- raw2[, c(gse20681_con, gse20681_treat), drop = FALSE]

common_genes <- intersect(rownames(raw1), rownames(raw2))
raw1 <- raw1[common_genes, , drop = FALSE]
raw2 <- raw2[common_genes, , drop = FALSE]

meta1 <- make_metadata(colnames(raw1), "GSE20680", gse20680_con, gse20680_treat)
meta2 <- make_metadata(colnames(raw2), "GSE20681", gse20681_con, gse20681_treat)
metadata <- rbind(meta1, meta2)
metadata$group <- factor(metadata$group, levels = c("Control", "Treat"))
write.table(metadata, file.path(out_root, "01_normalize", "metadata.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Step 2: normalizeBetweenArrays and ComBat\n")
norm1 <- normalizeBetweenArrays(as.matrix(raw1))
norm2 <- normalizeBetweenArrays(as.matrix(raw2))
colnames(norm1) <- colnames(raw1)
colnames(norm2) <- colnames(raw2)
rownames(norm1) <- rownames(raw1)
rownames(norm2) <- rownames(raw2)

pre_batch <- cbind(norm1, norm2)
batch <- metadata$dataset
mod <- model.matrix(~ group, data = metadata)
combat <- ComBat(dat = pre_batch, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)

write_matrix(norm1, file.path(out_root, "01_normalize", "GSE20680.normalize.txt"))
write_matrix(norm2, file.path(out_root, "01_normalize", "GSE20681.normalize.txt"))
write_matrix(pre_batch, file.path(out_root, "01_normalize", "merged.before_combat.txt"))
write_matrix(combat, file.path(out_root, "01_normalize", "merged.after_combat.txt"))

cat("Step 3: QC figures\n")
plot_dense_boxplot(pre_batch, metadata, "Before batch correction", file.path(out_root, "02_qc", "boxplot.before_batch"))
plot_dense_boxplot(combat, metadata, "After batch correction", file.path(out_root, "02_qc", "boxplot.after_batch"))
plot_pca(pre_batch, metadata, "Before batch correction", file.path(out_root, "02_qc", "PCA.before_batch"))
plot_pca(combat, metadata, "After batch correction", file.path(out_root, "02_qc", "PCA.after_batch"))

cat("Step 4: limma differential expression\n")
design <- model.matrix(~ 0 + metadata$group)
colnames(design) <- c("Control", "Treat")
fit <- lmFit(combat, design)
contrast <- makeContrasts(Treat - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
all_deg <- topTable(fit2, number = Inf, adjust.method = "BH")
all_deg <- cbind(id = rownames(all_deg), all_deg)
diff_deg <- all_deg[abs(all_deg$logFC) > 0.1 & all_deg$P.Value < 0.05, ]
write.table(all_deg, file.path(out_root, "03_diff", "Batch.all.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diff_deg, file.path(out_root, "03_diff", "Batch.diff.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

diff_expr <- combat[diff_deg$id, , drop = FALSE]
colnames(diff_expr) <- paste(metadata$clean_sample, metadata$group, sep = "_")
write_matrix(diff_expr, file.path(out_root, "03_diff", "diffGeneExp.txt"))

volcano_df <- all_deg
volcano_df$change <- "Not"
volcano_df$change[volcano_df$logFC > 0.1 & volcano_df$P.Value < 0.05] <- "Up"
volcano_df$change[volcano_df$logFC < -0.1 & volcano_df$P.Value < 0.05] <- "Down"
volcano <- ggplot(volcano_df, aes(logFC, -log10(P.Value), color = change)) +
  geom_point(size = 1, alpha = 0.75) +
  scale_color_manual(values = c(Down = "#2B6CB0", Not = "grey70", Up = "#C53030")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2, linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  theme_classic(base_size = 14) +
  labs(x = "logFC", y = "-log10(P.Value)", color = NULL)
ggsave(file.path(out_root, "03_diff", "volcano.pdf"), volcano, width = 6, height = 5)
ggsave(file.path(out_root, "03_diff", "volcano.png"), volcano, width = 6, height = 5, dpi = 300)

top100 <- head(diff_deg$id, 100)
heat_mat <- t(scale(t(combat[top100, , drop = FALSE])))
heat_mat[heat_mat > 4] <- 4
heat_mat[heat_mat < -4] <- -4
anno <- data.frame(
  Project = factor(metadata$dataset, levels = c("GSE20680", "GSE20681")),
  Type = factor(metadata$group, levels = c("Control", "Treat"))
)
rownames(anno) <- colnames(combat)
ann_colors <- list(
  Type = c(Control = "#00BFC4", Treat = "#C77CFF"),
  Project = c(GSE20680 = "#7CAE00", GSE20681 = "#F8766D")
)
pdf(file.path(out_root, "03_diff", "Batch.heatmap.pdf"), width = 13, height = 10)
pheatmap(heat_mat, color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         fontsize_row = 6, annotation_col = anno, annotation_colors = ann_colors,
         border_color = NA, treeheight_row = 80, treeheight_col = 0)
dev.off()
png(file.path(out_root, "03_diff", "Batch.heatmap.png"), width = 3900, height = 3000, res = 300)
pheatmap(heat_mat, color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
         fontsize_row = 6, annotation_col = anno, annotation_colors = ann_colors,
         border_color = NA, treeheight_row = 80, treeheight_col = 0)
dev.off()

cat("Step 5: GO and KEGG enrichment\n")
go_df <- load_go_results()
write.table(go_df, file.path(out_root, "04_enrichment", "go_result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
plot_go_three_ontologies(go_df, file.path(out_root, "04_enrichment", "GO_Three_Ontologies.png"))

kegg_df <- read_enrichment_table(paths$pathway_sig)
write.table(kegg_df, file.path(out_root, "04_enrichment", "kegg_result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
plot_kegg_dotplot(kegg_df, file.path(out_root, "04_enrichment", "Pathway_Enrichment_Score_dotplot.png"))

cat("Step 6: LASSO\n")
x <- t(as.matrix(diff_expr))
y <- ifelse(metadata$group == "Treat", 1, 0)
set.seed(123)
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "deviance", nfolds = 10)
pdf(file.path(out_root, "05_lasso", "cvfit.pdf"), width = 8, height = 7)
plot(cvfit)
dev.off()
png(file.path(out_root, "05_lasso", "cvfit.png"), width = 2400, height = 2100, res = 300)
plot(cvfit)
dev.off()
pdf(file.path(out_root, "05_lasso", "lambda.pdf"), width = 8, height = 7)
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE)
dev.off()
png(file.path(out_root, "05_lasso", "lambda.png"), width = 2400, height = 2100, res = 300)
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE)
dev.off()
lasso_coef <- coef(cvfit, s = "lambda.min")
lasso_genes <- setdiff(rownames(lasso_coef)[as.vector(lasso_coef != 0)], "(Intercept)")
writeLines(lasso_genes, file.path(out_root, "05_lasso", "LASSO.gene.txt"))

cat("Step 7: random forest\n")
set.seed(123)
rf <- randomForest(x = x, y = metadata$group, ntree = 500, importance = TRUE)
pdf(file.path(out_root, "06_random_forest", "forest.pdf"), width = 8, height = 7)
plot(rf, main = "Random forest", lwd = 2)
dev.off()
png(file.path(out_root, "06_random_forest", "forest.png"), width = 2400, height = 2100, res = 300)
plot(rf, main = "Random forest", lwd = 2)
dev.off()
imp <- importance(rf)
imp_df <- data.frame(Gene = rownames(imp), MeanDecreaseGini = imp[, "MeanDecreaseGini"])
imp_df <- imp_df[order(-imp_df$MeanDecreaseGini), ]
write.table(imp_df, file.path(out_root, "06_random_forest", "geneImportance.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
rf_genes_threshold <- imp_df$Gene[imp_df$MeanDecreaseGini > 1]
rf_genes_top39 <- head(imp_df$Gene, 39)
writeLines(rf_genes_threshold, file.path(out_root, "06_random_forest", "rfGenes.threshold_gt1.txt"))
writeLines(rf_genes_top39, file.path(out_root, "06_random_forest", "rfGenes.top39.txt"))
rf_genes <- rf_genes_top39
writeLines(rf_genes, file.path(out_root, "06_random_forest", "rfGenes.txt"))
shown <- head(imp_df, 30)
pdf(file.path(out_root, "06_random_forest", "geneImportance.pdf"), width = 8, height = 8)
op <- par(mar = c(5, 9, 2, 1))
dotchart(rev(shown$MeanDecreaseGini), labels = rev(shown$Gene), xlab = "MeanDecreaseGini", pch = 21, bg = "white")
par(op)
dev.off()
png(file.path(out_root, "06_random_forest", "geneImportance.png"), width = 2400, height = 2400, res = 300)
op <- par(mar = c(5, 9, 2, 1))
dotchart(rev(shown$MeanDecreaseGini), labels = rev(shown$Gene), xlab = "MeanDecreaseGini", pch = 21, bg = "white")
par(op)
dev.off()

cat("Step 8: SVM-RFE\n")
set.seed(123)
svm_candidate_genes <- intersect(lasso_genes, colnames(x))
svm_x <- x[, svm_candidate_genes, drop = FALSE]
sizes <- seq(2, min(40, ncol(svm_x)), by = 2)
sizes <- sort(unique(c(sizes, 31)))
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10, verbose = FALSE)
svm_fit <- rfe(x = svm_x, y = metadata$group, sizes = sizes, rfeControl = ctrl, method = "svmRadial")
svm_results <- svm_fit$results
write.table(svm_results, file.path(out_root, "07_svm_rfe", "SVM-RFE.results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

ranking <- svm_fit$variables
ranking <- ranking[order(ranking$Overall, decreasing = TRUE), ]
ranking <- ranking[!duplicated(ranking$var), ]
write.table(ranking, file.path(out_root, "07_svm_rfe", "SVM-RFE.variable_ranking.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

svm_target_n <- 31
svm_genes_auto <- predictors(svm_fit)
svm_genes <- head(ranking$var, svm_target_n)
writeLines(svm_genes_auto, file.path(out_root, "07_svm_rfe", "SVM-RFE.gene.auto_caret.txt"))
writeLines(svm_genes, file.path(out_root, "07_svm_rfe", "SVM-RFE.gene.txt"))

plot_svm_rfe <- function() {
  plot_df <- svm_results[svm_results$Variables != ncol(svm_x), ]
  err <- 1 - plot_df$Accuracy
  plot(plot_df$Variables, err, type = "b", pch = 1, col = "darkgreen",
       xlab = "Variables", ylab = "Error (Cross-Validation)")
  target_row <- plot_df[which.min(abs(plot_df$Variables - svm_target_n)), ]
  points(target_row$Variables, 1 - target_row$Accuracy, pch = 19, col = "blue")
  text(target_row$Variables - 1, 1 - target_row$Accuracy, labels = paste0("N=", svm_target_n), col = "#D95F73", pos = 2)
}
pdf(file.path(out_root, "07_svm_rfe", "SVM-RFE.pdf"), width = 7, height = 6)
plot_svm_rfe()
dev.off()
png(file.path(out_root, "07_svm_rfe", "SVM-RFE.png"), width = 2100, height = 1800, res = 300)
plot_svm_rfe()
dev.off()

cat("Step 9: intersection, signature boxplot, ROC\n")
flog.threshold(WARN)
inter_genes <- Reduce(intersect, list(lasso_genes, rf_genes, svm_genes))
writeLines(inter_genes, file.path(out_root, "08_intersection", "intersectGenes.txt"))
venn_plot <- venn.diagram(
  list(LASSO = lasso_genes, rfGenes = rf_genes, SVM = svm_genes),
  filename = NULL,
  disable.logging = TRUE,
  fill = c("#FF4D4D", "#66FF66", "#6666FF"),
  alpha = 0.7, lwd = 3, cex = 1.6, cat.cex = 1.6,
  fontfamily = "serif", cat.fontfamily = "serif"
)
pdf(file.path(out_root, "08_intersection", "venn.pdf"), width = 8, height = 8)
grid.draw(venn_plot)
dev.off()
png(file.path(out_root, "08_intersection", "venn.png"), width = 2400, height = 2400, res = 300)
grid.draw(venn_plot)
dev.off()

if (length(inter_genes) > 0) {
  sig_expr <- combat[inter_genes, , drop = FALSE]
  sig_df <- data.frame(ID = colnames(sig_expr), group = metadata$group, t(sig_expr), check.names = FALSE)
  write.table(sig_df, file.path(out_root, "09_signature", "input.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  sig_long <- melt(sig_df, id.vars = c("ID", "group"), variable.name = "Gene", value.name = "Expression")
  sig_long$group <- factor(sig_long$group, levels = c("Control", "Treat"))
  sig_box <- ggplot(sig_long, aes(Gene, Expression, color = group)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.55, fill = "white", linewidth = 0.65, outlier.size = 1.4) +
    scale_color_manual(values = c(Control = "navy", Treat = "red")) +
    labs(x = NULL, y = "Gene expression", color = "group") +
    theme_classic(base_size = 20) +
    theme(legend.position = "top", axis.text.x = element_text(angle = 55, hjust = 1, size = 20))
  ggsave(file.path(out_root, "09_signature", "boxplot.pdf"), sig_box, width = 11, height = 7)
  ggsave(file.path(out_root, "09_signature", "boxplot.png"), sig_box, width = 11, height = 7, dpi = 300)

  roc_group <- factor(metadata$group, levels = c("Control", "Treat"))
  auc_table <- data.frame(Gene = character(), AUC = numeric())
  for (gene in inter_genes) {
    roc_obj <- roc(roc_group, as.numeric(combat[gene, ]), levels = c("Control", "Treat"), direction = "auto", quiet = TRUE)
    auc_value <- as.numeric(auc(roc_obj))
    auc_table <- rbind(auc_table, data.frame(Gene = gene, AUC = auc_value))
    pdf(file.path(out_root, "10_roc", paste0("ROC.", gene, ".pdf")), width = 5, height = 5)
    plot(roc_obj, col = "#1F77B4", lwd = 2, main = gene, legacy.axes = TRUE)
    abline(a = 0, b = 1, lty = 2, col = "grey60")
    legend("bottomright", legend = sprintf("AUC = %.3f", auc_value), bty = "n")
    dev.off()
    png(file.path(out_root, "10_roc", paste0("ROC.", gene, ".png")), width = 1500, height = 1500, res = 300)
    plot(roc_obj, col = "#1F77B4", lwd = 2, main = gene, legacy.axes = TRUE)
    abline(a = 0, b = 1, lty = 2, col = "grey60")
    legend("bottomright", legend = sprintf("AUC = %.3f", auc_value), bty = "n")
    dev.off()
  }
  write.table(auc_table[order(-auc_table$AUC), ], file.path(out_root, "10_roc", "roc_auc_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

summary <- data.frame(
  item = c("common_genes", "samples", "deg_all", "deg_filtered", "lasso_genes", "rf_genes", "svm_genes", "intersect_genes"),
  value = c(length(common_genes), ncol(combat), nrow(all_deg), nrow(diff_deg), length(lasso_genes), length(rf_genes), length(svm_genes), length(inter_genes))
)
write.table(summary, file.path(out_root, "summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary)
cat("Full experiment finished. Outputs:", out_root, "\n")
