# T-MDSP reproducibility repository

This directory contains only the material requested for public release:

1. Simulation 1 and Simulation 2 missing-40% twenty-repeat workflows.
2. Lung and brain T-MDSP training code.
3. Unsupervised PCA-dimension selection (Horn parallel analysis).
4. PCA-only five-cluster evaluation with ARI and NMI.
5. Small reference tables used by `scripts/verify_reproduction.R`.

The lung and brain analyses use exactly the same penalty definition:

`kappa_ij(lambda) = kappa0 * ((1-lambda) + lambda*w_ij)`.

Only the configuration changes: lung uses `lambda=1`, so it is the fully
adaptive-Z model; brain uses `lambda=0.5`. The coefficient is therefore not a
different algorithm in the two datasets.

Dataset-specific starting states are stored under `data/initialization/` and
are loaded by the corresponding training configuration.

## Run

From this directory:

```r
source("scripts/install_packages.R")
source("scripts/run_all.R")
```

For a quick code/data audit without the long model fits:

```r
source("scripts/verify_reproduction.R")
```

Raw lung and brain count matrices are included under `data/real/`. Simulation
inputs are generated locally, so the multi-hundred-MB intermediate RDS files
from the working analysis are deliberately excluded from GitHub.
