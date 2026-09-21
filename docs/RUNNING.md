# Running the analysis

## Inspect the included results

Run `Rscript scripts/check_results.R` from the root. Only base R and the published result tables are used. The script confirms that the three saved gene lists yield the saved intersection, checks that the external AUC table covers the same genes, and prints its values. It does not validate the scientific methods or rerun experiments.

## Inputs

The 28 project inputs are included through Git LFS. Install Git LFS and run `git lfs pull` after cloning if necessary. Run `Rscript code/01_check_environment.R` to verify the files, dimensions, gene lists, and package availability.

## Install and run

```sh
Rscript code/install_packages.R
Rscript code/run_final.R
```

The installer obtains the CRAN and Bioconductor dependencies listed in `code/00_config.R`. Historical package versions were not recorded, and no historical lockfile is supplied. The saved-input evaluation has been run successfully with R 4.6.0 and the current dependency set.

To run the separate experimental path, use `Rscript code/run_full_experiment.R`. It writes to `generated/20_full_experiment/` and does not overwrite the saved-input evaluation. The GO and KEGG steps plot the included enrichment tables; they do not reconstruct the historical enrichment query from scratch.

| Script | Purpose |
|---|---|
| `00_config.R` | Paths, package lists, and shared helpers |
| `01_check_environment.R` | Check local files and report package availability |
| `20_full_pipeline.R` | Separate experimental path for normalization, differential analysis and feature screening |
| `30_generate_figures.R` | Generate figures using saved intermediate inputs, with some model refitting |
| `40_external_validation.R` | Probe mapping, external single-gene ROC and combined scores |
| `run_final.R` | Run the saved-input evaluation, figure generation, external evaluation, and run metadata |
| `run_full_experiment.R` | Run the separate normalization, differential-expression, and feature-selection path |

The separate experimental path does not automatically feed its outputs into `run_final.R`. Both paths read existing enrichment tables. Figure generation also reads saved feature lists; refitting a model for a plot does not regenerate every published list.

## Outputs

New outputs go to `generated/`. Scripts clear their corresponding generated subfolders before writing, so keep manually created files elsewhere. Published snapshots in `results/` and `figures/` remain separate.

## Verification status

Repository preparation included R syntax checks, missing-input handling, and checks of the included result tables and documentation links. The saved-input evaluation has been rerun after restoring the inputs. The corrected external-validation code freezes training-cohort direction and scaling; its generated summary is separate from the historical result snapshot. The fresh experimental path also runs, but its current-environment DEG and intersection counts differ from the historical snapshot; see `generated/20_full_experiment/RECOMPUTATION_NOTE.txt` after running it. See [changes](CHANGELOG.md).
