rm(list = ls())
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese_China.936"))

# Final Simulation 1 missing40 repeated-mask generator.
#
# This script encodes the final mask mechanism selected for the B=20 repeat:
# 1. fixed rho_true / sim_cov / labels / site_type from Simulation 1;
# 2. replay the original simulation generator to recover unsaved intermediates
#    such as sim_cov_mean, latent_score, marker_owner, base_missing, and
#    missing_score;
# 3. for each replicate, regenerate missing20 by the original Bernoulli +
#    adjust-to-target logic;
# 4. generate nested missing40 from that replicate-specific missing20 mask by
#    the original missing40 score rule;
# 5. calibrate endpoint-marker missing rate to 0.685 by coverage-bin-conserving
#    swaps into background / continuous_lowrank entries, without restoring any
#    missing20 entry.
#
# Recommended Windows/R invocation:
#   cd "C:/Users/86186/Desktop/tensor/dna data/结果/simulation1_f"
#   Rscript repeat_missing40_B20/code/make_final_missing40_masks_B20.R

set.seed(20260704)

if (!file.exists(file.path("main", "inputs", "missing20_source.rds"))) {
  stop(
    "Run this script from the simulation1_f root directory, e.g. ",
    "Rscript repeat_missing40_B20/code/make_final_missing40_masks_B20.R"
  )
}
root_dir <- "."
input_dir <- file.path(root_dir, "main", "inputs")
out_dir <- file.path(root_dir, "repeat_missing40_B20", "mask_generation_final")
mask_dir <- file.path(out_dir, "replicate_masks")
dir.create(mask_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("SIM_REPEAT_B", unset = "20"))
G <- 5L
k_rank <- 5L

input20 <- readRDS(file.path(input_dir, "missing20_source.rds"))
input40 <- readRDS(file.path(input_dir, "missing40_source.rds"))

rho_true <- input20$rho_true
sim_cov <- input20$sim_cov
labels <- input20$labels
site_type <- input20$site_type
original_mask20 <- input20$mask
original_mask40 <- input40$mask

stopifnot(identical(rho_true, input40$rho_true))
stopifnot(identical(sim_cov, input40$sim_cov))
stopifnot(identical(labels, input40$labels))
stopifnot(identical(site_type, input40$site_type))

n_sites <- nrow(rho_true)
n_cells <- ncol(rho_true)
n_entries <- n_sites * n_cells
target_missing20 <- round(0.20 * n_entries)
target_missing40 <- round(0.40 * n_entries)

sigmoid <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))
normalize_rows <- function(X, target = 1) X / sqrt(rowSums(X^2) + 1e-8) * target

match_exact01_by_class <- function(rho, labels, target_zero, target_one) {
  G <- length(target_zero)
  rho_out <- rho
  for (g in seq_len(G)) {
    cells_g <- which(labels == g)
    idx_g <- which(col(rho) %in% cells_g)
    n_g <- length(idx_g)
    n_zero <- max(0L, min(n_g, round(n_g * target_zero[g])))
    n_one <- max(0L, min(n_g - n_zero, round(n_g * target_one[g])))
    ord <- idx_g[order(rho[idx_g], decreasing = FALSE)]
    if (n_zero > 0) {
      rho_out[ord[seq_len(n_zero)]] <- 0
    }
    if (n_one > 0) {
      one_idx <- ord[seq.int(n_g - n_one + 1L, n_g)]
      rho_out[one_idx] <- 1
    }
  }
  rho_out
}

adjust_mask_to_target <- function(mask, score, target_missing) {
  delta <- target_missing - sum(!mask)
  if (delta > 0) {
    observed_idx <- which(mask)
    add_idx <- observed_idx[order(score[observed_idx], decreasing = TRUE)[seq_len(delta)]]
    mask[add_idx] <- FALSE
  } else if (delta < 0) {
    missing_idx <- which(!mask)
    restore_idx <- missing_idx[order(score[missing_idx], decreasing = FALSE)[seq_len(-delta)]]
    mask[restore_idx] <- TRUE
  }
  matrix(mask, nrow(mask), ncol(mask))
}

replay_original_components <- function(setting = c("missing20", "missing40")) {
  setting <- match.arg(setting)
  set.seed(20260618)

  background_fraction <- 0.80
  endpoint_fraction <- 0.045
  continuous_fraction <- max(0.01, 1 - background_fraction - endpoint_fraction)
  endpoint_missing_rate <- 0.44
  endpoint_owner_missing_rate <- 0.44
  endpoint_nonowner_missing_rate <- 0.44
  continuous_missing_rate <- 0.15
  background_missing_rate <- 0.12
  marker_shift <- 2.0
  continuous_uv_weight <- 0.56
  endpoint_uv_weight <- 0.05
  v_noise_sd <- 0.96
  v_class_scale <- 1.64
  background_class_sd <- 0.22
  fuzzy_pair_enabled <- FALSE
  nuisance_enabled <- FALSE
  nonlinear_nuisance_enabled <- FALSE
  hierarchy_enabled <- FALSE

  labels_replay <- sample(rep(seq_len(G), length.out = n_cells))
  lung500_top5_prop <- c(168, 103, 81, 77, 39) / sum(c(168, 103, 81, 77, 39))
  class_sizes <- as.integer(floor(n_cells * lung500_top5_prop))
  class_sizes[seq_len(n_cells - sum(class_sizes))] <-
    class_sizes[seq_len(n_cells - sum(class_sizes))] + 1L
  labels_replay <- sample(rep(seq_len(G), times = class_sizes))

  R_true <- matrix(0, n_cells, G)
  R_true[cbind(seq_len(n_cells), labels_replay)] <- 1

  site_type_replay <- sample(
    c("background", "endpoint_marker", "continuous_lowrank"),
    n_sites,
    replace = TRUE,
    prob = c(background_fraction, endpoint_fraction, continuous_fraction)
  )

  M_true <- matrix(0, G, k_rank)
  for (g in seq_len(min(G, k_rank))) {
    M_true[g, g] <- v_class_scale
  }
  M_true <- scale(M_true, center = TRUE, scale = FALSE)
  V_true <- R_true %*% M_true + matrix(rnorm(n_cells * k_rank, sd = v_noise_sd), n_cells, k_rank)
  U_true <- normalize_rows(matrix(rnorm(n_sites * k_rank), n_sites, k_rank), sqrt(k_rank))

  bg_comp_prob <- c(0.09, 0.586, 0.16, 0.164)
  bg_comp <- sample(1:4, n_sites, replace = TRUE, prob = bg_comp_prob / sum(bg_comp_prob))
  base_eta <- numeric(n_sites)
  base_eta[bg_comp == 1] <- rnorm(sum(bg_comp == 1), -5.5, 0.35)
  base_eta[bg_comp == 2] <- rnorm(sum(bg_comp == 2), 5.5, 0.35)
  base_eta[bg_comp == 3] <- rnorm(sum(bg_comp == 3), 0, 0.35)
  base_eta[bg_comp == 4] <- rnorm(sum(bg_comp == 4), 0, 1.00)

  eta_class <- matrix(base_eta, n_sites, G)
  marker_owner <- sample(seq_len(G), n_sites, replace = TRUE)
  marker_direction <- sample(c(-1, 1), n_sites, replace = TRUE)
  for (g in seq_len(G)) {
    idx <- site_type_replay == "endpoint_marker" & marker_owner == g
    eta_class[idx, g] <- marker_shift * marker_direction[idx]
    eta_class[idx, setdiff(seq_len(G), g)] <- -marker_shift * marker_direction[idx]
  }

  background_class_effect <- matrix(rnorm(n_sites * G, sd = background_class_sd), n_sites, G)
  background_class_effect <- background_class_effect - matrix(rowMeans(background_class_effect), n_sites, G)
  eta_class[site_type_replay == "background", ] <- eta_class[site_type_replay == "background", ] +
    background_class_effect[site_type_replay == "background", ]

  uv_weight <- rep(0.04, n_sites)
  uv_weight[site_type_replay == "continuous_lowrank"] <- continuous_uv_weight
  uv_weight[site_type_replay == "endpoint_marker"] <- endpoint_uv_weight
  eta_true <- eta_class[, labels_replay] + (U_true %*% t(V_true)) * uv_weight +
    matrix(rnorm(n_sites * n_cells, sd = 0.22), n_sites, n_cells)

  lung500_global_methylation <- 0.791
  lung500_class_methylation <- c(0.805, 0.841, 0.687, 0.772, 0.812)
  lung500_class_shift_strength <- 0.65
  class_eta_shift <- lung500_class_shift_strength *
    (logit(lung500_class_methylation) - logit(lung500_global_methylation))
  eta_true <- eta_true + matrix(class_eta_shift[labels_replay], n_sites, n_cells, byrow = TRUE)

  global_eta_shift <- 0
  for (i in seq_len(12)) {
    current_mean <- mean(sigmoid(eta_true + global_eta_shift))
    step <- logit(lung500_global_methylation) - logit(current_mean)
    global_eta_shift <- global_eta_shift + step
    if (abs(step) < 1e-4) break
  }
  eta_true <- eta_true + global_eta_shift

  rho_replay <- sigmoid(eta_true)
  rho_replay <- match_exact01_by_class(
    rho_replay,
    labels_replay,
    target_zero = c(0.079840, 0.066162, 0.145291, 0.103380, 0.077044),
    target_one = c(0.598402, 0.666231, 0.433082, 0.564034, 0.609551)
  )

  base_missing <- matrix(background_missing_rate, n_sites, n_cells)
  base_missing[site_type_replay == "endpoint_marker", ] <- endpoint_missing_rate
  base_missing[site_type_replay == "continuous_lowrank", ] <- continuous_missing_rate
  if (any(site_type_replay == "endpoint_marker")) {
    endpoint_sites <- which(site_type_replay == "endpoint_marker")
    for (g in seq_len(G)) {
      owner_site_g <- endpoint_sites[marker_owner[endpoint_sites] == g]
      cells_g <- which(labels_replay == g)
      other_cells <- setdiff(seq_len(n_cells), cells_g)
      if (length(owner_site_g) > 0 && length(cells_g) > 0) {
        base_missing[owner_site_g, cells_g] <- endpoint_owner_missing_rate
      }
      if (length(owner_site_g) > 0 && length(other_cells) > 0) {
        base_missing[owner_site_g, other_cells] <- endpoint_nonowner_missing_rate
      }
    }
  }

  latent_score <- scale(U_true %*% t(V_true))
  class_missing_multiplier <- c(1.18, 0.92, 1.08, 1.00, 1.28)
  base_missing <- base_missing * matrix(class_missing_multiplier[labels_replay], n_sites, n_cells, byrow = TRUE)

  site_cov <- rlnorm(n_sites, meanlog = log(3.3), sdlog = 0.45)
  cell_cov <- rlnorm(n_cells, meanlog = log(4.2), sdlog = 0.22)
  sim_cov_mean <- outer(site_cov, cell_cov, "*") / mean(cell_cov)
  sim_cov_replay <- matrix(rpois(n_sites * n_cells, lambda = pmax(sim_cov_mean, 0.05)), n_sites, n_cells)
  coverage_low_score <- as.numeric(scale(-log1p(sim_cov_mean)))
  coverage_low_score <- matrix(coverage_low_score, n_sites, n_cells)

  base_missing <- base_missing + 0.10 * coverage_low_score
  extreme01 <- rho_replay == 0 | rho_replay == 1
  base_missing <- base_missing + 0.20 * extreme01
  base_missing <- pmin(pmax(base_missing + 0.08 * as.numeric(latent_score > 0.8), 0.02), 0.85)

  missing_score <- base_missing + 0.05 * coverage_low_score
  if (setting == "missing40") {
    missing_score <- missing_score + 0.12 * matrix(site_type_replay == "continuous_lowrank", n_sites, n_cells)
    missing_score <- missing_score + 0.14 * extreme01
    missing_score <- missing_score - 0.03 * matrix(site_type_replay == "background", n_sites, n_cells)
  }

  list(
    rho_true = rho_replay,
    labels = labels_replay,
    site_type = site_type_replay,
    sim_cov = sim_cov_replay,
    sim_cov_mean = sim_cov_mean,
    marker_owner = marker_owner,
    latent_score = latent_score,
    base_missing = base_missing,
    missing_score = missing_score
  )
}

original20 <- replay_original_components("missing20")
original40 <- replay_original_components("missing40")

stopifnot(isTRUE(all.equal(original20$rho_true, rho_true)))
stopifnot(identical(as.integer(original20$labels), as.integer(labels)))
stopifnot(identical(as.character(original20$site_type), as.character(site_type)))
stopifnot(identical(original20$sim_cov, sim_cov))
stopifnot(isTRUE(all.equal(original40$rho_true, rho_true)))
stopifnot(identical(original40$sim_cov, sim_cov))

make_coverage_bins <- function(coverage) {
  cov_values <- as.vector(coverage)
  cov_rank <- rank(cov_values, ties.method = "average")
  matrix(pmin(3L, pmax(1L, ceiling(3 * cov_rank / length(cov_rank)))), nrow(coverage), ncol(coverage))
}

endpoint_rate <- function(mask) {
  endpoint_mat <- matrix(site_type == "endpoint_marker", n_sites, n_cells)
  mean((!mask)[endpoint_mat])
}

calibrate_endpoint_missing <- function(mask40, mask20, target_high = 0.69, target_mid = 0.685,
                                       rep_id = 1) {
  current_rate <- endpoint_rate(mask40)
  if (current_rate <= target_high) {
    return(mask40)
  }

  set.seed(51000 + rep_id)
  mask <- mask40
  cov_bin <- make_coverage_bins(sim_cov)
  endpoint_mat <- matrix(site_type == "endpoint_marker", n_sites, n_cells)
  non_endpoint_mat <- !endpoint_mat
  site_mat <- matrix(site_type, n_sites, n_cells)
  n_endpoint_entries <- sum(endpoint_mat)
  n_restore <- ceiling((current_rate - target_mid) * n_endpoint_entries)
  restore_pool <- which(endpoint_mat & (!mask) & mask20)
  if (!length(restore_pool) || n_restore <= 0) {
    return(mask)
  }
  n_restore <- min(n_restore, length(restore_pool))

  restore_bins <- cov_bin[restore_pool]
  bin_counts <- table(factor(restore_bins, levels = 1:3))
  alloc <- floor(n_restore * as.numeric(bin_counts) / sum(bin_counts))
  remainder <- n_restore - sum(alloc)
  if (remainder > 0) {
    frac <- n_restore * as.numeric(bin_counts) / sum(bin_counts) - alloc
    ord <- order(frac, decreasing = TRUE)
    alloc[ord[seq_len(remainder)]] <- alloc[ord[seq_len(remainder)]] + 1L
  }
  alloc <- pmin(alloc, as.numeric(bin_counts))

  restore_idx <- integer(0)
  for (g in seq_len(3)) {
    pool_g <- restore_pool[restore_bins == g]
    if (alloc[g] > 0 && length(pool_g) > 0) {
      restore_idx <- c(restore_idx, sample(pool_g, alloc[g]))
    }
  }
  if (!length(restore_idx)) {
    return(mask)
  }
  mask[restore_idx] <- TRUE

  restored_bins <- cov_bin[restore_idx]
  original_site_rates <- tapply((!original_mask40), matrix(site_type, n_sites, n_cells), mean)
  site_entries <- table(factor(site_type, levels = c("background", "continuous_lowrank", "endpoint_marker")))
  target_bg_count <- round(original_site_rates[["background"]] * site_entries[["background"]])
  target_cont_count <- round(original_site_rates[["continuous_lowrank"]] * site_entries[["continuous_lowrank"]])
  current_bg_count <- sum((!mask) & (site_mat == "background"))
  current_cont_count <- sum((!mask) & (site_mat == "continuous_lowrank"))
  deficit <- c(
    background = max(0, target_bg_count - current_bg_count),
    continuous_lowrank = max(0, target_cont_count - current_cont_count)
  )
  if (sum(deficit) == 0) {
    deficit[] <- c(1, 1)
  }

  add_idx <- integer(0)
  for (g in seq_len(3)) {
    need_g <- sum(restored_bins == g)
    if (need_g <= 0) next

    bg_pool <- which((site_mat == "background") & (cov_bin == g) & mask)
    cont_pool <- which((site_mat == "continuous_lowrank") & (cov_bin == g) & mask)
    n_bg <- min(length(bg_pool), round(need_g * deficit[["background"]] / sum(deficit)))
    n_cont <- need_g - n_bg
    if (n_cont > length(cont_pool)) {
      n_bg <- min(length(bg_pool), n_bg + n_cont - length(cont_pool))
      n_cont <- min(n_cont, length(cont_pool))
    }
    if (n_bg > length(bg_pool)) {
      n_cont <- min(length(cont_pool), n_cont + n_bg - length(bg_pool))
      n_bg <- min(n_bg, length(bg_pool))
    }

    chosen <- integer(0)
    if (n_bg > 0) chosen <- c(chosen, sample(bg_pool, n_bg))
    if (n_cont > 0) chosen <- c(chosen, sample(cont_pool, n_cont))
    if (length(chosen) < need_g) {
      fallback_pool <- setdiff(which(non_endpoint_mat & (cov_bin == g) & mask), chosen)
      n_extra <- min(length(fallback_pool), need_g - length(chosen))
      if (n_extra > 0) chosen <- c(chosen, sample(fallback_pool, n_extra))
    }
    if (length(chosen) < need_g) {
      fallback_pool <- setdiff(which(non_endpoint_mat & mask), chosen)
      n_extra <- min(length(fallback_pool), need_g - length(chosen))
      if (n_extra > 0) chosen <- c(chosen, sample(fallback_pool, n_extra))
    }

    mask[chosen] <- FALSE
    add_idx <- c(add_idx, chosen)
  }

  stopifnot(sum(!mask) == sum(!mask40))
  stopifnot(!any((!mask20) & mask))
  mask
}

mask_site_diagnostics <- function(mask, rep_id) {
  do.call(rbind, lapply(sort(unique(site_type)), function(tp) {
    idx <- matrix(site_type == tp, n_sites, n_cells)
    data.frame(
      replicate = rep_id,
      total_missing_rate = mean(!mask),
      site_type = tp,
      site_type_missing_rate = mean((!mask)[idx]),
      site_type_missing_fraction = sum((!mask)[idx]) / sum(!mask),
      stringsAsFactors = FALSE
    )
  }))
}

mask_coverage_diagnostics <- function(mask, rep_id) {
  cov_bin <- make_coverage_bins(sim_cov)
  data.frame(
    replicate = rep_id,
    coverage_group = c("low", "middle", "high"),
    missing_rate = vapply(seq_len(3), function(g) mean((!mask)[cov_bin == g]), numeric(1)),
    entries = vapply(seq_len(3), function(g) sum(cov_bin == g), numeric(1)),
    stringsAsFactors = FALSE
  )
}

mask_label_diagnostics <- function(mask, rep_id) {
  data.frame(
    replicate = rep_id,
    true_label = sort(unique(labels)),
    missing_rate = vapply(sort(unique(labels)), function(lb) {
      mean((!mask)[, labels == lb, drop = FALSE])
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

make_final_masks <- function(rep_id) {
  set.seed(41000 + rep_id)
  mask20_initial <- matrix(runif(n_entries) > original20$base_missing, n_sites, n_cells)
  score20 <- original20$missing_score * runif(n_entries)
  mask20 <- adjust_mask_to_target(mask20_initial, score20, target_missing20)

  set.seed(42000 + rep_id)
  invisible(runif(n_entries) > original40$base_missing)
  score40 <- original40$missing_score * runif(n_entries)
  mask40_raw <- adjust_mask_to_target(mask20, score40, target_missing40)
  mask40 <- calibrate_endpoint_missing(mask40_raw, mask20, rep_id = rep_id)

  stopifnot(sum(!mask20) == target_missing20)
  stopifnot(sum(!mask40) == target_missing40)
  stopifnot(all(!mask40[!mask20]))

  list(mask20 = mask20, mask40_raw = mask40_raw, mask40 = mask40)
}

all_site <- list()
all_cov <- list()
all_label <- list()
manifest <- vector("list", B)

for (b in seq_len(B)) {
  cat(sprintf("[%s] make final missing40 mask %d/%d\n", format(Sys.time(), "%H:%M:%S"), b, B))
  masks <- make_final_masks(b)
  saveRDS(
    list(
      replicate = b,
      mask20 = masks$mask20,
      mask40_raw_before_endpoint_calibration = masks$mask40_raw,
      mask40 = masks$mask40,
      rho_true = rho_true,
      sim_cov = sim_cov,
      labels = labels,
      site_type = site_type,
      protocol = "original replay missing20 + nested original missing40 + endpoint rate 0.685 coverage-bin-conserving calibration"
    ),
    file.path(mask_dir, sprintf("replicate_%02d_masks.rds", b))
  )
  all_site[[b]] <- mask_site_diagnostics(masks$mask40, b)
  all_cov[[b]] <- mask_coverage_diagnostics(masks$mask40, b)
  all_label[[b]] <- mask_label_diagnostics(masks$mask40, b)
  manifest[[b]] <- data.frame(
    replicate = b,
    missing20_rate = mean(!masks$mask20),
    missing40_raw_rate = mean(!masks$mask40_raw),
    missing40_rate = mean(!masks$mask40),
    endpoint_raw_rate = endpoint_rate(masks$mask40_raw),
    endpoint_final_rate = endpoint_rate(masks$mask40),
    nested = all(!masks$mask40[!masks$mask20]),
    stringsAsFactors = FALSE
  )
}

site_diag <- do.call(rbind, all_site)
cov_diag <- do.call(rbind, all_cov)
label_diag <- do.call(rbind, all_label)
manifest <- do.call(rbind, manifest)

write.csv(manifest, file.path(out_dir, "final_mask_manifest.csv"), row.names = FALSE)
write.csv(site_diag, file.path(out_dir, "final_mask_site_type_diagnostics.csv"), row.names = FALSE)
write.csv(cov_diag, file.path(out_dir, "final_mask_coverage_diagnostics.csv"), row.names = FALSE)
write.csv(label_diag, file.path(out_dir, "final_mask_label_diagnostics.csv"), row.names = FALSE)

summary_site <- aggregate(cbind(site_type_missing_rate, site_type_missing_fraction) ~ site_type, site_diag, mean)
summary_cov <- aggregate(missing_rate ~ coverage_group, cov_diag, mean)
summary_label <- aggregate(missing_rate ~ true_label, label_diag, mean)

write.csv(summary_site, file.path(out_dir, "final_mask_site_type_summary.csv"), row.names = FALSE)
write.csv(summary_cov, file.path(out_dir, "final_mask_coverage_summary.csv"), row.names = FALSE)
write.csv(summary_label, file.path(out_dir, "final_mask_label_summary.csv"), row.names = FALSE)

readme <- c(
  "# Final Simulation 1 missing40 B=20 mask generator",
  "",
  "This folder is generated by `repeat_missing40_B20/code/make_final_missing40_masks_B20.R`.",
  "",
  "Final mechanism:",
  "1. Fixed `rho_true`, `sim_cov`, `labels`, and `site_type` from `main/inputs/missing20_source.rds`.",
  "2. Replayed the original Simulation 1 generator with `set.seed(20260618)` and the original missing20/missing40 parameters to recover unsaved `base_missing`, `missing_score`, `sim_cov_mean`, `marker_owner`, and `latent_score`.",
  "3. For each replicate b, regenerated missing20 using `set.seed(41000 + b)`, the original Bernoulli mask, and the original adjust-to-20% rule.",
  "4. Regenerated nested missing40 using `set.seed(42000 + b)`, the original missing40 score rule, and no restoration of missing20 entries.",
  "5. If endpoint-marker missing rate exceeded 0.69, restored endpoint-marker entries added after missing20 until the endpoint rate reached 0.685, then added the same number of missing values to background/continuous entries within the same coverage bin.",
  "",
  "Primary output:",
  "- `replicate_masks/replicate_XX_masks.rds`: contains `mask20`, raw nested `mask40`, final calibrated `mask40`, and fixed data identifiers.",
  "- `final_mask_manifest.csv`: per-replicate missing rates and nested checks.",
  "- `final_mask_*_diagnostics.csv`: site-type, coverage-bin, and label diagnostics.",
  "- `final_mask_*_summary.csv`: B-replicate averages."
)
writeLines(readme, file.path(out_dir, "README_final_mask_generation.md"))

cat("\nSite-type summary:\n")
print(summary_site)
cat("\nCoverage summary:\n")
print(summary_cov)
cat("\nManifest endpoint rates:\n")
print(manifest[, c("replicate", "missing40_rate", "endpoint_raw_rate", "endpoint_final_rate", "nested")])
