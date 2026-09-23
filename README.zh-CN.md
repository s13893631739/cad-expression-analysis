# CAD 转录组分析与预测签名项目

这是一个可独立展示、可复现运行的冠心病（CAD）基因表达分析项目，使用 R 完成从数据校验到预测签名评估的完整流程。

## 项目概览

项目整合两个发现队列和一个外部队列，完成质量控制、标准化、批次校正、差异表达、通路富集、多方法候选筛选、嵌套 Elastic Net 建模和外部评估。

## 主要结果

- **发现数据：** GSE20680 和 GSE20681 共 337 个样本、19,261 个共同基因。
- **差异表达：** `abs(logFC) > 0.1` 且名义 `P < 0.05` 得到 286 个基因；其中 1 个通过 BH FDR < 0.05。
- **候选筛选：** LASSO 32 个、随机森林 39 个、SVM-RFE 31 个，描述性交集 15 个。
- **预测模型：** 嵌套 Elastic Net 外层 OOF ROC-AUC **0.5505**、PR-AUC **0.5969**，最终模型 6 个非零基因。
- **外部评估：** GSE113079 冻结模型探索性评估 ROC-AUC **0.4415**、PR-AUC **0.6481**。

主要文件：

- [项目结果摘要](results/project_summary.md)
- [指标表](results/project_summary.tsv)
- [预测签名系数](results/predictive_signature.tsv)
- [方法说明](docs/METHODS.md)
- [运行说明](docs/RUNNING.md)

项目图表：

- [批次校正后 PCA](figures/pca_after_batch.png)
- [差异表达火山图](figures/volcano.png)
- [GO 富集](figures/go_enrichment.png)
- [嵌套模型 ROC](figures/nested_model_roc.png)
- [外部评估 ROC](figures/external_model_roc.png)

## 数据

| 数据集 | 用途 | CAD | 对照 |
|---|---|---:|---:|
| GSE20680 | 发现队列 | 87 | 52 |
| GSE20681 | 发现队列 | 99 | 99 |
| GSE113079 | 外部评估 | 93 | 48 |

原始 GEO 矩阵和样本分组文件通过 Git LFS 管理；外部探针注释使用 GPL20115。

## 分析流程

1. 校验输入文件和 SHA-256；
2. 标准化发现队列、合并共同基因并使用 ComBat 校正数据集批次；
3. 生成 QC 图并进行 limma CAD/Control 差异分析；
4. 使用明确的基因背景和多重校正重算 GO/KEGG；
5. 进行 LASSO、随机森林和 SVM-RFE 描述性候选筛选；
6. 使用嵌套 Elastic Net 建模，在训练折内完成特征筛选、中心化/标准化和 alpha 选择；
7. 固定最终模型系数，在 GSE113079 上进行外部评估。

## 一键运行

```sh
git lfs install
git lfs pull
Rscript scripts/check_inputs.R
Rscript code/run_project.R
```

完整输出写入 `generated/project/`，可审阅的核心结果提交在 `results/`。

## 项目局限

这是探索性研究流程，不是经过临床验证的诊断工具。外部队列用于探索性评估，不是盲法确认队列；跨平台表达可比性、概率校准、置信区间和完全未查看的确认队列仍需进一步验证。

## 简历描述

使用两个 GEO 队列的 337 个样本构建可复现 CAD 转录组 R 分析流程，完成标准化、ComBat 批次校正、limma 差异分析、GO/KEGG 富集和多方法候选筛选；在训练折内完成特征筛选与标准化，使用嵌套 Elastic Net 得到 OOF ROC-AUC 0.5505、PR-AUC 0.5969 的 6 基因预测模型，并打包输入校验、报告和一键运行入口。

[English README](README.md)
