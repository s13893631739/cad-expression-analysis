# Running the analysis

Run the commands below from the repository root. The first checks are read-only and use base R; the analysis steps may install packages and write files under `generated/`.

## 1. Clone and download the data

The repository stores 28 large inputs (8 raw files and 20 intermediate files) with Git LFS. Install Git LFS once on the machine, then run:

```sh
git lfs install
git lfs pull
```

Check that the expected objects are present:

```sh
git lfs ls-files
Rscript scripts/check_inputs.R
```

`check_inputs.R` reads `data/input_manifest.tsv` and checks every listed file for existence, a hydrated (non-pointer) LFS object, byte size, and SHA-256 checksum. It does not train a model.

If this check reports a missing file or an LFS pointer, run `git lfs pull` again. If you cloned without LFS installed, install it and run `git lfs pull` from the repository root.

## 2. Inspect the saved result tables

```sh
Rscript scripts/check_results.R
```

This base-R check confirms that the saved LASSO, random-forest, and SVM-RFE lists yield the saved intersection, and that the saved external AUC table covers the same genes. A successful run reports 32, 39, and 31 genes in the three lists, 10 shared candidates, and the saved AUC values. It does not retrain models or recompute ROC curves.

## 3. Check the R environment

```sh
Rscript code/01_check_environment.R
```

This reports the R version, data and output paths, basic input dimensions, candidate-list counts, and whether all required CRAN and Bioconductor packages are installed. It stops if a required input is missing or is still an LFS pointer.

The scripts support alternate locations when needed:

```sh
CAD_DATA_ROOT=/path/to/data CAD_OUTPUT_ROOT=/path/to/output \
  Rscript code/01_check_environment.R
```

`CAD_DATA_ROOT` must contain the same `raw/` and `intermediate/` layout as the repository's `data/` directory.

## 4. Install packages and run the saved-input evaluation

```sh
Rscript code/install_packages.R
Rscript code/run_final.R
```

The installer obtains the CRAN and Bioconductor dependencies listed in `code/00_config.R`. Historical package versions were not recorded, and no lockfile is supplied. The saved-input evaluation has been run successfully with R 4.6.0 and the current dependency set.

`run_final.R` evaluates the restored intermediate inputs, generates figures and external-validation summaries, and records `sessionInfo()` and run metadata. New files are written below `generated/`; the historical snapshots in `results/` and `figures/` are kept separate.

## 5. Run the separate experimental path (optional)

```sh
Rscript code/run_full_experiment.R
```

This runs the normalization, differential-expression, and feature-selection path from the restored inputs. It writes to `generated/20_full_experiment/` and does not overwrite the saved-input evaluation. The current-environment counts can differ from the historical snapshot; read `generated/20_full_experiment/RECOMPUTATION_NOTE.txt` after the run.

The GO and KEGG steps plot the included enrichment tables; they do not reconstruct the historical enrichment query from GEO. The separate experimental path does not automatically feed its outputs into `run_final.R`. Figure generation also reads saved feature lists; refitting a model for a plot does not regenerate every published list.

## Output locations

- `results/`: published numerical snapshots and gene lists kept unchanged.
- `figures/`: selected historical figures.
- `generated/`: outputs from the current run, including `RUN_METADATA.txt` and `sessionInfo.txt`.
- `generated/20_full_experiment/`: outputs and the comparison note for the optional fresh experimental path.

Scripts clear their corresponding generated subfolders before writing, so keep manually created files elsewhere.

## Reproduction scope

The saved-input pipeline is runnable from the restored files. A complete GEO-download-to-input reconstruction is not included: the original acquisition, annotation, grouping, and preprocessing steps used to create some intermediate inputs were not fully recorded. Current reruns therefore provide a reproducible execution of the supplied inputs, while historical result snapshots remain explicitly separate from newly computed outputs. See [the change log](CHANGELOG.md) and [the data notes](../data/README.md).
