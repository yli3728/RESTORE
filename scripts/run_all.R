scripts <- c(
  "simulation/simulation1/generate_data_and_masks.R",
  "simulation/simulation1/run_B20.R",
  "simulation/simulation2/generate_data_and_masks.R",
  "simulation/simulation2/run_B20.R",
  "real_data/train_lung.R",
  "real_data/select_pca_and_cluster.R --dataset=lung",
  "real_data/train_brain.R",
  "real_data/select_pca_and_cluster.R --dataset=brain",
  "scripts/verify_reproduction.R"
)
message("Run these commands from the repository root:")
for (x in scripts) message("Rscript ", x)
