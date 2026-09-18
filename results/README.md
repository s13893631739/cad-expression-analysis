# Saved results

These files are copied from the original project and have not been recomputed during repository preparation.

| File | Contents |
|---|---|
| [LASSO list](05_lasso/LASSO.gene.txt) | 32 candidate genes |
| [Random forest list](06_random_forest/rfGenes.txt) | 39 candidate genes |
| [SVM-RFE list](07_svm/SVM-RFE.gene.txt) | 31 candidate genes |
| [SVM-RFE results](07_svm/SVM-RFE.results.txt) | Saved internal cross-validation metrics |
| [Intersection](08_venn/intersectGenes.txt) | 10 shared candidates |
| [Internal AUC](10_roc/roc_auc_table.txt) | Single-gene discovery-set AUC |
| [External AUC](11_external_validation/roc_auc_table.txt) | Single-gene AUC in GSE113079 |
| [Probe mapping](11_external_validation/probe_mapping.txt) | External probe-to-gene mapping |

## Figures

- [Candidate overlap](../figures/feature_selection.png)
- [External single-gene ROC](../figures/external_single_gene_roc.png)
- [Volcano plot](../figures/volcano.png), using the unadjusted P-value threshold
- [GO enrichment](../figures/go_enrichment.png), based on existing enrichment tables

The historical weighted-score AUC of 0.539 is transcribed from the original `external_validation_summary.txt`. That summary also contained the affected transferred-GLM metric, so it is not reproduced as a current model evaluation table. See [methods](../docs/METHODS.md) for interpretation.
