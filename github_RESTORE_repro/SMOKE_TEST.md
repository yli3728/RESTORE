# Verification record

The release was audited on 2026-09-14 with R 4.5.2 on Windows.

- Every public R script parses successfully.
- Both real-data trainers completed a one-iteration end-to-end run.
- The clustering script completed for both datasets with a two-permutation PA smoke setting.
- The formula audit confirmed a single `kappa_lambda()` implementation.
- Lung resolved to `lambda=1`; brain resolved to `lambda=0.5`.

The complete 200-iteration training and 100-permutation PCA selection were run
for both datasets. After
freezing each dataset's complete continuous initialization state, the public
workflow reproduced the reference conclusions:

| dataset | selected PCs | new ARI | new NMI | reference ARI | reference NMI |
|---|---:|---:|---:|---:|---:|
| lung | 6 | 0.9333570 | 0.9052022 | 0.9333570 | 0.9052022 |
| brain | 5 | 0.9703429 | 0.9428095 | 0.9703429 | 0.9428095 |

The two datasets load their own complete starting states.

The legacy targets are retained in
`reference_results/lung_clustering_reference.csv`,
`reference_results/brain_clustering_reference.csv`, and
`reference_results/simulation_summary.csv`.
