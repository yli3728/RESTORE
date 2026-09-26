rm(list = ls())
set.seed(2026)

required <- c("Matrix", "RSpectra", "aricode")
for (p in required) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)

sigmoid <- function(x) 1 / (1 + exp(-x))
normalize_rows <- function(X, target = 1) X / sqrt(rowSums(X^2) + 1e-8) * target
simplex_row <- function(v) {
  u <- sort(v, decreasing = TRUE); cssv <- cumsum(u) - 1; ind <- seq_along(u)
  rho <- max(ind[u - cssv / ind > 0]); pmax(v - cssv[rho] / rho, 0)
}
project_simplex <- function(M) t(apply(M, 1, simplex_row))
mean_fill <- function(Y, mask) {
  P <- Y; mu <- rowSums(Y * mask) / (rowSums(mask) + 1e-8)
  for (i in seq_len(nrow(P))) P[i, !mask[i, ]] <- mu[i]
  P
}

run_lambda050 <- function(Y_obs, mask, coverage, mode, k_rank = 5, G = 5,
                         max_iter_override = NULL) {
  # All Simulation 1 lambda-Z experiments use this fixed coefficient.
  z_lambda <- 0.5
  n_sites <- nrow(Y_obs); n_cells <- ncol(Y_obs)
  if (mode == "simulation1") {
    max_iter <- 180L; main_prob <- 0.85; kappa <- 1.0
    lambda_3 <- 0.8; lr_U <- 0.0012; lr_M <- 0.0012; lr_R <- 0.0002
    lambda_l2 <- 0.01; lambda_s <- 0.04; trunc_start <- 70L
    lower <- 0.04; upper <- 0.96; dual_lr <- 0.12; normalize_uv <- TRUE
  } else {
    max_iter <- 200L; main_prob <- 0.78; kappa <- 0.8
    lambda_3 <- 1.2; lr_U <- 0.001; lr_M <- 0.001; lr_R <- 0.0001
    lambda_l2 <- 0.01; lambda_s <- 0.04; trunc_start <- 100L
    lower <- 0.02; upper <- 0.95; dual_lr <- 0.1; normalize_uv <- FALSE
  }
  if (!is.null(max_iter_override)) max_iter <- as.integer(max_iter_override)
  W <- Matrix::bandSparse(n_sites, k = 1, diagonals = list(rep(1, n_sites - 1)), symmetric = TRUE)
  Dinv <- Matrix::Diagonal(x = 1 / (Matrix::rowSums(W) + 1e-8))
  L <- Matrix::Diagonal(n_sites) - Dinv %*% W
  rm(W, Dinv)

  Ybase <- mean_fill(Y_obs, mask)
  sv <- RSpectra::svds(scale(t(Ybase), center = TRUE, scale = FALSE), k = k_rank)
  Vinit <- sv$u %*% diag(sv$d, nrow = k_rank)
  U <- normalize_rows(sv$v, sqrt(k_rank))
  set.seed(2026); km <- kmeans(Vinit, centers = G, nstart = 60)
  M <- normalize_rows(km$centers, sqrt(k_rank))
  R <- matrix((1 - main_prob) / (G - 1), n_cells, G)
  R[cbind(seq_len(n_cells), km$cluster)] <- main_prob
  site_init <- rowSums(Y_obs * mask) / (rowSums(mask) + 1e-8)
  s <- qlogis(pmin(pmax(site_init, 0.02), 0.98))
  Z <- if (mode == "simulation1") Ybase else matrix(0.5, n_sites, n_cells)
  Lambda <- matrix(0, n_sites, n_cells)
  rm(Ybase, sv, Vinit, km); gc()

  # The entry-specific coverage proxy uses smooth site/cell summaries,
  # including at entries hidden by the mask; its mixing coefficient is fixed above.
  site_cov <- rowMeans(coverage); cell_cov <- colMeans(coverage)
  site_scale <- median(site_cov[site_cov > 0]); cell_scale <- median(cell_cov[cell_cov > 0])
  coverage_factor <- outer(site_cov / (site_scale + 1e-8),
                           cell_cov / (cell_scale + 1e-8), "+")
  coverage_factor <- coverage_factor / (coverage_factor + 1)
  support <- pmin(1, rowSums(mask) / (median(rowSums(mask)[rowSums(mask) > 0]) + 1e-8))

  start <- Sys.time()
  for (iter in seq_len(max_iter)) {
    V <- R %*% M
    eta <- pmin(pmax(U %*% t(V) + matrix(s, n_sites, n_cells), -12), 12)
    rho <- sigmoid(eta)
    polar <- pmax(0, 1 - 4 * rho * (1 - rho))
    entry_weight <- polar * coverage_factor * support
    kappa_weight <- (1 - z_lambda) + z_lambda * entry_weight
    grad <- (rho - Y_obs) * mask + kappa * kappa_weight * (rho - Z + Lambda)
    U <- U - lr_U * pmin(pmax(grad %*% V + lambda_l2 * U, -5), 5)
    U <- scale(U, center = TRUE, scale = FALSE)
    if (normalize_uv) U <- normalize_rows(U, sqrt(k_rank))
    gradV <- t(grad) %*% U
    M <- M - lr_M * pmin(pmax(t(R) %*% gradV + lambda_l2 * M, -5), 5)
    M <- scale(M, center = TRUE, scale = FALSE)
    if (normalize_uv) M <- normalize_rows(M, sqrt(k_rank))
    R <- project_simplex(R - lr_R * pmin(pmax(gradV %*% t(M), -5), 5))
    grads <- rowSums(grad) + lambda_3 * as.numeric(L %*% s) + lambda_s * s
    grads <- grads / (rowSums(mask) + 1e-8)
    s <- pmin(pmax(s - 0.001 * pmin(pmax(grads, -5), 5), -8), 8)
    X <- rho + Lambda
    Z <- if (iter > trunc_start) ifelse(X < lower, 0, ifelse(X > upper, 1, X)) else rho
    Z <- pmin(pmax(Z, 0), 1)
    Lambda <- pmin(pmax(Lambda + dual_lr * (rho - Z), -0.5), 0.5)
    if (iter %% 20 == 0) cat(sprintf("  lambda=0.5 Z iter %d/%d, elapsed %.1fs\n", iter, max_iter,
                                     as.numeric(Sys.time() - start, units = "secs")))
  }
  V <- R %*% M
  fixed <- U %*% t(V)
  shift <- rowSums(fixed * mask) / (rowSums(mask) + 1e-8)
  fixed <- fixed - matrix(shift, n_sites, n_cells); s <- pmin(pmax(s + shift, -8), 8)
  for (ii in 1:25) {
    rho <- sigmoid(pmin(pmax(fixed + matrix(s, n_sites, n_cells), -12), 12))
    step <- rowSums((rho - Y_obs) * mask) / (rowSums(mask * rho * (1 - rho)) + 1e-8)
    step <- pmin(pmax(step, -1), 1); s <- pmin(pmax(s - step, -8), 8)
    if (mean(abs(step)) < 1e-5) break
  }
  P <- sigmoid(pmin(pmax(fixed + matrix(s, n_sites, n_cells), -12), 12))
  list(P_hat = P, V_hat = V, coverage_weight_scale = 1, z_lambda = z_lambda,
       max_iter = max_iter)
}

