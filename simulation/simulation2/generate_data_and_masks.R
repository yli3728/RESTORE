rm(list = ls())

if (!file.exists(file.path("main", "inputs", "restored_all_missing20_40_60_results.rds"))) {
  stop("Run from the simulation2_f root directory")
}

input_path <- file.path("main", "inputs", "restored_all_missing20_40_60_results.rds")
out_root <- file.path("repeat_missing40_B20", "glm_site_gumbel_9500_masks_B20")
mask_dir <- file.path(out_root, "replicate_masks")
dir.create(mask_dir, recursive = TRUE, showWarnings = FALSE)

shared <- readRDS(input_path)
rho_true <- shared$rho_true
coverage <- shared$cov_mtx
truth <- shared$true_labels
site_type <- shared$site_type
mask20 <- shared$results$missing20$mask
reference40 <- shared$results$missing40$mask
n_sites <- nrow(rho_true)
n_cells <- ncol(rho_true)
n_total <- length(rho_true)

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(.5, length(x)))
  (x - rng[1]) / diff(rng)
}
clip <- function(x) pmin(pmax(x, 1e-8), 1 - 1e-8)

coverage_cell_score <- 1 - rescale01(colMeans(coverage))
coverage_entry_score <- matrix(coverage_cell_score, n_sites, n_cells, byrow = TRUE)
continuous_entry <- matrix(site_type == "continuous_lowrank", n_sites, n_cells)
background_entry <- matrix(site_type == "background", n_sites, n_cells)
endpoint_entry <- matrix(site_type == "endpoint_marker", n_sites, n_cells)
extreme01_entry <- rho_true <= .01 | rho_true >= .99
reference_add <- mask20 & !reference40

# Reproduce the frozen candidate fit exactly.
set.seed(20261000L)
fit_idx <- sample.int(n_total, min(750000L, n_total))
fit_data <- data.frame(
  added = as.integer(reference_add[fit_idx]),
  cov_low = as.vector(coverage_entry_score)[fit_idx],
  endpoint = as.integer(as.vector(endpoint_entry)[fit_idx]),
  continuous = as.integer(as.vector(continuous_entry)[fit_idx]),
  background = as.integer(as.vector(background_entry)[fit_idx]),
  extreme01 = as.integer(as.vector(extreme01_entry)[fit_idx]),
  class = factor(truth[((fit_idx - 1L) %/% n_sites) + 1L])
)
fit_added <- glm(added ~ cov_low + endpoint + continuous + background + extreme01 + class,
                 data = fit_data, family = binomial())
prediction_data <- data.frame(
  cov_low = as.vector(coverage_entry_score),
  endpoint = as.integer(as.vector(endpoint_entry)),
  continuous = as.integer(as.vector(continuous_entry)),
  background = as.integer(as.vector(background_entry)),
  extreme01 = as.integer(as.vector(extreme01_entry)),
  class = factor(rep(truth, each = n_sites), levels = levels(fit_data$class))
)
lp <- predict(fit_added, newdata = prediction_data, type = "link")
add_probability <- matrix(clip(plogis(lp)), n_sites, n_cells)

additional_by_site <- pmax(rowSums(!reference40) - rowSums(!mask20), 0L)
if (sum(additional_by_site) != sum(!reference40) - sum(!mask20)) {
  stop("Per-site additional counts do not sum to the missing40 extension target")
}

make_mask <- function(rep_id) {
  set.seed(9500L + rep_id)
  score <- log(add_probability) - log(-log(runif(n_total)))
  mask <- mask20
  for (i in seq_len(n_sites)) {
    add <- additional_by_site[i]
    if (add <= 0L) next
    candidates <- which(mask[i, ])
    chosen <- candidates[order(score[i, candidates], decreasing = TRUE)[seq_len(add)]]
    mask[i, chosen] <- FALSE
  }
  mask
}

manifest <- vector("list", 20L)
for (b in seq_len(20L)) {
  cat(sprintf("Write glm-site-gumbel mask %d/20\n", b))
  mask40 <- make_mask(b)
  stopifnot(all(!mask40[!mask20]), sum(!mask40) == sum(!reference40),
            identical(rowSums(!mask40), rowSums(!reference40)))
  mask_path <- file.path(mask_dir, sprintf("replicate_%02d_missing40_glm_site_gumbel.rds", b))
  saveRDS(list(
    simulation = "simulation2", replicate = b, seed = 9500L + b,
    mask20 = mask20, mask40 = mask40,
    protocol = paste0("fixed missing20; fitted add probability from coverage/site_type/",
                      "extreme01/class; exact reference additional count per site; ",
                      "probability-weighted Gumbel ranking")
  ), mask_path)
  manifest[[b]] <- data.frame(
    replicate = b, seed = 9500L + b, missing_count = sum(!mask40),
    missing_rate = mean(!mask40), missing20_preserved = all(!mask40[!mask20]),
    exact_reference_per_site_counts = identical(rowSums(!mask40), rowSums(!reference40)),
    mask_file = basename(mask_path), mask_md5 = unname(tools::md5sum(mask_path))
  )
  rm(mask40); gc(verbose = FALSE)
}
manifest <- do.call(rbind, manifest)
write.csv(manifest, file.path(out_root, "mask_manifest.csv"), row.names = FALSE)
saveRDS(list(
  input_md5 = unname(tools::md5sum(input_path)), fit_seed = 20261000L,
  seed0 = 9500L, coefficients = coef(fit_added),
  additional_by_site = additional_by_site, add_probability_summary = summary(as.vector(add_probability))
), file.path(out_root, "protocol_payload.rds"))

writeLines(c(
  "# Simulation2 GLM site-Gumbel missing40 masks (B=20)", "",
  "Approximate Simulation1-style recovery candidate selected by MeanFill validation.",
  "It preserves the saved missing20 base and exact reference missing40 counts per site,",
  "then randomizes the selected cells using fitted probabilities and Gumbel ranking.", "",
  "Historical MeanFill recovery is approximate, not exact."
), file.path(out_root, "README.md"))
print(manifest)
