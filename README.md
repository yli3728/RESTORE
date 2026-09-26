# RESTORE: Structure-preserving reconstruction of sparse single-cell DNA methylomes

<p align="center">
  <img src="restore-logo.png" alt="RESTORE logo" width="360">
</p>

RESTORE reconstructs sparse single-cell DNA methylomes while retaining the biological structure needed for downstream cell-population analysis. It models methylation on the probability scale with a low-rank cell representation, region-specific effects, graph smoothing across neighboring genomic regions, and an adaptive auxiliary-matrix penalty that becomes stronger for confident methylation states.

This repository contains the R implementation and reproducibility materials for RESTORE: two 40%-missing simulation studies, lung and brain single-cell methylation analyses, fixed initialization states, reference results, and verification scripts. For the real datasets, RESTORE reconstructs a region-by-cell methylation-probability matrix and then evaluates cell structure using Horn parallel analysis, PCA, k-means clustering, adjusted Rand index (ARI), and normalized mutual information (NMI).

## Workflow

1. Prepare matched region-by-cell methylated-read and coverage matrices.
2. Convert observed counts to methylation proportions and retain a coverage mask.
3. Fit RESTORE to obtain the reconstructed methylation matrix `P_hat` and cell representation `V_hat`.
4. Select the number of principal components by 95th-percentile permutation Horn parallel analysis.
5. Cluster cells with k-means and, when reference labels are available, report ARI and NMI.
6. Verify the model formula, dataset settings, initialization objects, and reproduced results against the checked-in references.

## Installation

RESTORE requires R and the packages `Matrix`, `aricode`, and `irlba`. The simulation workflows additionally use `RSpectra`.

From an R session opened in the repository root, install the core dependencies with:

```r
source("scripts/install_packages.R")
install.packages("RSpectra") # required for the simulation workflows
```

Alternatively:

```r
install.packages(c("Matrix", "aricode", "irlba", "RSpectra"))
```

No package build step is required; the analysis is run directly from the repository root.

## Input data

### Real-data input

Each real dataset is represented by two RDS files with identical dimensions and dimnames:

| Object | Required format |
|---|---|
| Coverage | Numeric region-by-cell matrix; zero denotes no observed coverage |
| Methylated reads | Numeric region-by-cell matrix of methylated-read counts |
| Rows | Genomic regions/sites, with stable row names |
| Columns | Cells, with stable cell IDs |
| Missing values | Accepted and converted to zero before the observation mask is constructed |

The observed methylation proportion is `methy / cov` wherever `cov > 0`. Dataset-specific initialization files in `data/initialization/` must match the input cell IDs, site IDs, and configured rank.

The included datasets are:

| Dataset | Coverage | Methylated reads | Cells | 
|---|---|---|---:|
| Lung | `data/real/LG_ACCPU_5kb_573cell_cov.rds` | `data/real/LG_ACCPU_5kb_573cell_methy.rds` | 573 | 
| Brain | `data/real/M1C_H1930001_5kb_542cell_cov.rds` | `data/real/M1C_H1930001_5kb_542cell_methy.rds` | 542 | 
### Simulation input

The simulation workflows use a true methylation-probability matrix, coverage information, cell labels, and replicate-specific masks with exactly 40% missing entries. The observed matrix is passed to RESTORE; the truth at held-out entries is used only for reconstruction evaluation.

## Quick start

Run a short, non-fitting audit from the repository root:

```r
source("scripts/verify_reproduction.R")
```

To fit and evaluate one included real dataset from a terminal:

```bash
Rscript real_data/train_lung.R
Rscript real_data/select_pca_and_cluster.R --dataset=lung
```

For the brain dataset, replace `lung` with `brain`:

```bash
Rscript real_data/train_brain.R
Rscript real_data/select_pca_and_cluster.R --dataset=brain
```

The fitted model is written to `outputs/<dataset>/model.rds`; clustering metrics and labels are written alongside it.

## Main functions

### `kappa_lambda()`

Defined in `R/model_formula.R`, this is the single implementation of RESTORE's adaptive penalty:

```r
kappa_lambda(kappa0, lambda, adaptive_weight)
```

`kappa0` is the base penalty, `lambda` controls the fixed-to-adaptive interpolation, and `adaptive_weight` is a non-negative region-by-cell weight.

### `train_tmdsp()`

Defined in `real_data/train_common.R`, this is the shared real-data fitting routine used by both dataset entry points:

```r
source("real_data/train_common.R")
train_tmdsp("lung")  # or "brain"
```

The routine loads the dataset configuration and fixed initialization state, constructs the observed methylation matrix and mask, optimizes the low-rank and structured components, and saves:

| Component | Meaning |
|---|---|
| `P_hat` | Reconstructed region-by-cell methylation probabilities |
| `V_hat` | Learned cell representation |
| `cell_ids` | Cell order associated with `P_hat` and `V_hat` |
| `lambda` | Dataset-specific adaptive-penalty setting |
| `formula` | Recorded penalty formula |

### `select_pca_and_cluster.R`

This script performs 100-permutation Horn parallel analysis by default, embeds cells using the selected PCs, runs five-cluster k-means, and writes ARI/NMI metrics plus cell-level assignments. Set `PA_B` to reduce the number of permutations for a smoke test:

```bash
PA_B=2 Rscript real_data/select_pca_and_cluster.R --dataset=lung
```

### Simulation routines

Each simulation directory contains a missingness generator, the RESTORE model core, and a 20-replicate runner. The runners evaluate held-out reconstruction with MSE, RMSE, and MAE and cell-structure recovery with ARI and NMI.

## Backward compatibility

The public method name is **RESTORE**, but the reproducibility code retains the earlier internal names `T-MDSP`, `train_tmdsp()`, `TMDSP_MAX_ITER`, `TMDSP_B20_OUTPUT_DIR`, and historical output filenames. These names are intentionally preserved so that archived configurations, reference tables, and existing reproduction commands continue to work. They refer to the RESTORE implementation in this repository.

## Included example data and output

The repository includes both input data and compact reference outputs:

| Path | Contents |
|---|---|
| `data/real/` | Lung and brain coverage and methylated-read matrices |
| `data/initialization/` | Dataset-specific model states and evaluation labels |
| `outputs/lung/` | Included lung labels, metrics, and training parameters |
| `outputs/brain/` | Included brain labels, metrics, and training parameters |
| `reference_results/` | Frozen real-data and simulation reference summaries |

The checked-in real-data results are:

| Dataset | Selected PCs | ARI | NMI |
|---|---:|---:|---:|
| Lung | 6 | 0.9333570 | 0.9052022 |
| Brain | 5 | 0.9703429 | 0.9428095 |

The simulation reference tables summarize 20 replicates at 40% missingness and include reconstruction error, clustering accuracy, variability, and comparisons with baseline methods.

## Reproduce the complete example

From a terminal in the repository root, run the real-data analyses in sequence:

```bash
Rscript real_data/train_lung.R
Rscript real_data/select_pca_and_cluster.R --dataset=lung
Rscript real_data/train_brain.R
Rscript real_data/select_pca_and_cluster.R --dataset=brain
Rscript scripts/verify_reproduction.R
```

The full fits use 200 optimization iterations and the clustering analyses use 100 parallel-analysis permutations, so they can take substantially longer than the audit. For a fitting smoke test, set `TMDSP_MAX_ITER=1`; for a clustering smoke test, set `PA_B=2`.

`scripts/run_all.R` prints the complete ordered command list, including both simulation workflows:

```r
source("scripts/run_all.R")
```

Simulation source objects and generated replicate masks are intentionally not duplicated in the compact release. To rerun those experiments, retain the directory layout expected by the scripts in `simulation/simulation1/` and `simulation/simulation2/`; the frozen 20-replicate summaries remain available in `reference_results/`.

## Repository structure

```text
R/
  model_formula.R                 Shared adaptive-penalty definition

config/
  datasets.csv                   Dataset-specific ranks, lambda values, and run settings

data/
  real/                          Lung and brain count matrices
  initialization/                Fixed model states and reference cell labels

real_data/
  train_common.R                 Shared RESTORE fitting implementation
  train_lung.R                   Lung entry point
  train_brain.R                  Brain entry point
  select_pca_and_cluster.R       Parallel analysis, PCA, clustering, ARI, and NMI

simulation/
  simulation1/                   Simulation 1 generator, model core, and B=20 runner
  simulation2/                   Simulation 2 generator, model core, and B=20 runner

scripts/
  install_packages.R             Core dependency installer
  run_all.R                      Ordered reproduction command list
  verify_reproduction.R          Formula, configuration, state, and result audit

outputs/
  lung/                          Included lung outputs
  brain/                         Included brain outputs

reference_results/               Frozen real-data and simulation reference tables
SMOKE_TEST.md                    Verification record and reproduced conclusions
restore-logo.png                RESTORE project logo
```
