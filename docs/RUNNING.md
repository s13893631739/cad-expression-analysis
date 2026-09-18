# Running the analysis

## Inspect the included results

Run `Rscript scripts/check_results.R` from the root. Only base R and the published result tables are used. The script confirms that the three saved gene lists yield the saved intersection, checks that the external AUC table covers the same genes, and prints its values. It does not validate the scientific methods or rerun experiments.

## Restore local inputs

Follow [data/README.md](../data/README.md). The environment check stops when files are missing. Obtaining the GEO records alone does not reconstruct all project-specific intermediate inputs.

## Install and run

```sh
Rscript code/install_packages.R
Rscript code/01_check_environment.R
Rscript code/run_final.R
```

The installer obtains the CRAN and Bioconductor dependencies listed in `code/00_config.R`. Historical package versions were not recorded, and no historical lockfile is supplied. Full execution with the current dependencies has not been verified.

| Script | Purpose |
|---|---|
| `00_config.R` | Paths, package lists, and shared helpers |
| `01_check_environment.R` | Check local files and report package availability |
| `20_full_pipeline.R` | Separate experimental path for normalization, differential analysis and feature screening |
| `30_generate_figures.R` | Generate figures using saved intermediate inputs, with some model refitting |
| `40_external_validation.R` | Probe mapping, external single-gene ROC and combined scores |
| `run_final.R` | Run environment check, figure generation, then external evaluation |

The separate experimental path does not automatically feed its outputs into `run_final.R`. Both paths read existing enrichment tables. Figure generation also reads saved feature lists; refitting a model for a plot does not regenerate every published list.

## Outputs

New outputs go to `generated/`. Scripts clear their corresponding generated subfolders before writing, so keep manually created files elsewhere. Published snapshots in `results/` and `figures/` remain separate.

## Verification status

Repository preparation included R syntax checks, missing-input handling, and checks of the included result tables and documentation links. The complete analysis was not rerun. The transferred GLM scoring correction therefore has no newly computed performance result. See [changes](CHANGELOG.md).
