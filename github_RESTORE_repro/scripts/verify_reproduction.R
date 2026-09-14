files <- c("R/model_formula.R",
  "simulation/simulation1/generate_data_and_masks.R",
  "simulation/simulation1/model_core.R", "simulation/simulation1/run_B20.R",
  "simulation/simulation2/generate_data_and_masks.R",
  "simulation/simulation2/model_core.R", "simulation/simulation2/run_B20.R",
  "real_data/train_lung.R",
  "real_data/train_brain.R", "real_data/select_pca_and_cluster.R")
for (f in files) parse(f)
cfg <- read.csv("config/datasets.csv")
stopifnot(cfg$lambda[cfg$dataset=="lung"] == 1,
          cfg$lambda[cfg$dataset=="brain"] == .5)
code <- paste(readLines("real_data/train_common.R"), collapse="\n")
stopifnot(grepl("kappa_lambda\\(kappa0,lambda,adaptive_weight\\)", code))
for (d in c("lung","brain")) {
  init <- readRDS(file.path("data/initialization", paste0(d, "_state.rds")))
  stopifnot(all(c("U_hat","s_hat","M_hat","R_mat") %in% names(init)))
  result <- file.path("outputs",d,"clustering_metrics.csv")
  if(file.exists(result)) {
    got <- read.csv(result); ref <- read.csv(file.path("reference_results",paste0(d,"_clustering_reference.csv")))
    stopifnot(abs(got$ARI-ref$ARI)<1e-7, abs(got$NMI-ref$NMI)<1e-7,
              got$PCA_dim==ref$PCA_dim)
  }
}
message("PASS: scripts parse; formula and lambdas are unified; dataset states are valid; available final results reproduce their references.")
