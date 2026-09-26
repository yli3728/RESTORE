rm(list = ls())

suite <- "simulation1"
B <- 20L
max_iter <- 180L
fit_seed0 <- 20260718L
evaluation_seed <- 2026L
evaluation_protocol <- paste0(
  "centered matrix PCA/SVD scores u*d; fixed_dim=3; kmeans K=5;",
  "nstart=50;iter.max=200;SVD_seed=2026;kmeans_seed=2026"
)

required <- c("Matrix", "RSpectra", "irlba", "aricode")
for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing R package: ", p)
}

# Run from this code directory. Keeping all R paths ASCII-only and relative
# avoids the Windows R locale corrupting the Chinese workspace path.
code_dir <- "."
repeat_dir <- ".."
root_dir <- file.path("..", "..")

mask_dir <- file.path(repeat_dir, "mask_generation_final", "replicate_masks")
default_out_dir <- file.path(repeat_dir, "results")
out_dir <- Sys.getenv("TMDSP_B20_OUTPUT_DIR", unset = default_out_dir)
if (!grepl("^([A-Za-z]:|/)", out_dir)) out_dir <- file.path(root_dir, out_dir)
dir.create(file.path(out_dir, "replicate_results"), recursive = TRUE, showWarnings = FALSE)

core_file <- file.path(code_dir, "model_core.R")
core_md5 <- unname(tools::md5sum(core_file))
ee <- parse(core_file)
needed_core <- c(
  "sigmoid", "normalize_rows", "simplex_row", "project_simplex",
  "mean_fill", "run_lambda050"
)
for (expr in ee) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.name(expr[[2]]) && as.character(expr[[2]]) %in% needed_core) {
    eval(expr, envir = .GlobalEnv)
  }
}
if (!all(vapply(needed_core, exists, logical(1), envir = .GlobalEnv, inherits = FALSE))) {
  stop("Could not load all T-MDSP core functions")
}

src <- readRDS(file.path(root_dir, "main", "inputs", "missing40_source.rds"))
rho_true <- src$rho_true
coverage <- src$sim_cov
truth <- as.integer(src$labels)
n_cells <- ncol(rho_true)

# Byte-for-byte equivalent to the frozen baseline PAC definition.
pac_dim3 <- function(P, labels) {
  X <- scale(t(P), center = TRUE, scale = FALSE)
  X[!is.finite(X)] <- 0
  set.seed(evaluation_seed)
  s <- irlba::irlba(X, nu = 3, nv = 3, maxit = 1000, work = 30)
  embedding <- sweep(s$u, 2, s$d, "*")
  set.seed(evaluation_seed)
  cluster <- kmeans(
    embedding, centers = 5, nstart = 50, iter.max = 200
  )$cluster
  list(
    ari = aricode::ARI(as.integer(cluster), as.integer(labels)),
    nmi = aricode::NMI(as.integer(cluster), as.integer(labels)),
    cluster = cluster,
    embedding = embedding
  )
}

rows <- vector("list", B)
for (bb in seq_len(B)) {
  cat(sprintf("\n[%s] lambdaZ=0.5, final mask40, PCA3 repeat %d/%d\n", suite, bb, B))
  mask_file <- file.path(mask_dir, sprintf("replicate_%02d_masks.rds", bb))
  if (!file.exists(mask_file)) stop("Missing final mask file: ", mask_file)
  mask_obj <- readRDS(mask_file)
  mask <- mask_obj$mask40
  stopifnot(
    identical(dim(mask), dim(rho_true)),
    identical(mask_obj$rho_true, rho_true),
    identical(mask_obj$sim_cov, coverage),
    identical(as.integer(mask_obj$labels), truth),
    sum(!mask) == round(0.40 * length(mask))
  )
  mask_md5 <- unname(tools::md5sum(mask_file))
  fit_seed <- fit_seed0 + bb
  set.seed(fit_seed)
  Y <- rho_true
  Y[!mask] <- 0
  fit <- run_lambda050(
    Y, mask, coverage, suite, max_iter_override = max_iter
  )
  pac <- pac_dim3(fit$P_hat, truth)
  idx <- !mask
  err <- fit$P_hat[idx] - rho_true[idx]
  rows[[bb]] <- data.frame(
    Simulation = suite, Replicate = bb, Missing_rate = mean(idx),
    Z_lambda = fit$z_lambda, Max_iter = fit$max_iter, Fit_seed = fit_seed,
    ARI = pac$ari, NMI = pac$nmi, Missing_MSE = mean(err^2),
    Missing_RMSE = sqrt(mean(err^2)), Missing_MAE = mean(abs(err)),
    Evaluation_protocol = evaluation_protocol,
    Mask_file = basename(mask_file), Mask_file_md5 = mask_md5,
    Core_code_md5 = core_md5, stringsAsFactors = FALSE
  )

  rep_dir <- file.path(out_dir, "replicate_results", sprintf("replicate_%02d", bb))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    data.frame(
      Cell = seq_len(n_cells), Truth = truth, Cluster = pac$cluster,
      PC1 = pac$embedding[, 1], PC2 = pac$embedding[, 2], PC3 = pac$embedding[, 3]
    ),
    file.path(rep_dir, "pca3_cluster_and_embedding.csv"), row.names = FALSE
  )
  saveRDS(
    list(
      P_hat = fit$P_hat, V_hat = fit$V_hat, z_lambda = fit$z_lambda,
      max_iter = fit$max_iter, fit_seed = fit_seed, cluster = pac$cluster,
      pca3_embedding = pac$embedding, truth = truth,
      evaluation_protocol = evaluation_protocol,
      mask_file = basename(mask_file), mask_file_md5 = mask_md5,
      core_code_md5 = core_md5, metrics = rows[[bb]]
    ),
    file.path(rep_dir, "lambda050Z_pca3_result.rds"), compress = TRUE
  )
  write.csv(
    do.call(rbind, rows[seq_len(bb)]),
    file.path(out_dir, "lambda050Z_pca3_B20_results_partial.csv"), row.names = FALSE
  )
  cat(sprintf(
    "  completed: ARI=%.6f, NMI=%.6f, missing MSE=%.8f, mask md5=%s\n",
    pac$ari, pac$nmi, rows[[bb]]$Missing_MSE, mask_md5
  ))
  rm(mask_obj, mask, Y, fit, pac, idx, err)
  gc()
}

res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "lambda050Z_pca3_B20_results.csv"), row.names = FALSE)
summary_table <- data.frame(
  Simulation = suite, B = B, Missing_rate = unique(res$Missing_rate),
  Z_lambda = unique(res$Z_lambda), Max_iter = unique(res$Max_iter),
  ARI_mean = mean(res$ARI), ARI_sd = sd(res$ARI),
  ARI_min = min(res$ARI), ARI_max = max(res$ARI),
  NMI_mean = mean(res$NMI), NMI_sd = sd(res$NMI),
  MSE_mean = mean(res$Missing_MSE), MSE_sd = sd(res$Missing_MSE),
  MSE_min = min(res$Missing_MSE), MSE_max = max(res$Missing_MSE),
  Evaluation_protocol = evaluation_protocol, Core_code_md5 = core_md5,
  stringsAsFactors = FALSE
)
write.csv(
  summary_table, file.path(out_dir, "lambda050Z_pca3_B20_summary.csv"), row.names = FALSE
)
writeLines(
  c(
    "# Simulation 1 T-MDSP missing40 B20", "",
    "- Masks: mask_generation_final/replicate_masks, field mask40.",
    "- Model: T-MDSP lambdaZ = 0.5, 180 iterations.",
    paste0("- Clustering: ", evaluation_protocol, "."),
    "- Each replicate RDS retains P_hat, V_hat, PCA3 embedding, labels, metrics, and file fingerprints."
  ),
  file.path(out_dir, "README.md")
)
print(summary_table)
