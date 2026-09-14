args <- commandArgs(trailingOnly=TRUE)
dataset <- sub("^--dataset=", "", args[grepl("^--dataset=", args)])
if (!dataset %in% c("lung", "brain")) stop("Use --dataset=lung or --dataset=brain")
if (!requireNamespace("aricode", quietly=TRUE)) stop("Install aricode")
b <- readRDS(file.path("outputs", dataset, "model.rds"))
lambda_value <- if (length(b$lambda)) b$lambda else if (dataset == "lung") 1 else .5
truth_tab <- read.csv(file.path("data/initialization", paste0(dataset, "_truth_labels.csv")))
truth <- factor(truth_tab$merged5[match(b$cell_ids, truth_tab$cell_id)])
keep <- !is.na(truth)
x <- scale(t(b$P_hat), center=TRUE, scale=FALSE); n <- nrow(x)
eig <- eigen(tcrossprod(x), symmetric=TRUE, only.values=FALSE)
observed <- pmax(eig$values, 0)
B <- as.integer(Sys.getenv("PA_B", "100")); null <- matrix(NA_real_, B, length(observed))
set.seed(2026)
for (bb in seq_len(B)) {
  set.seed(if (dataset == "lung") 920000L + bb else 930000L + bb)
  xp <- apply(x, 2, sample)
  null[bb, ] <- pmax(eigen(tcrossprod(xp), symmetric=TRUE, only.values=TRUE)$values, 0)
}
threshold <- apply(null, 2, quantile, .95)
p_dim <- max(1L, sum(cummax(observed <= threshold) == 0L))
scores <- sweep(eig$vectors[, seq_len(p_dim), drop=FALSE], 2, sqrt(observed[seq_len(p_dim)]), "*")
cluster_input <- if (dataset == "lung") scale(scores, center=TRUE, scale=TRUE) else scores
if (any(!is.finite(cluster_input))) stop("Non-finite PCA clustering input")
set.seed(2026); cluster <- kmeans(cluster_input, 5L, nstart=50L, iter.max=100L)$cluster
metrics <- data.frame(dataset=dataset, lambda=lambda_value, PCA_dim=p_dim,
  PCA_rule="Horn parallel analysis, permutation q95",
  ARI=aricode::ARI(as.integer(truth[keep]), cluster[keep]),
  NMI=aricode::NMI(as.integer(truth[keep]), cluster[keep]))
write.csv(metrics, file.path("outputs", dataset, "clustering_metrics.csv"), row.names=FALSE)
write.csv(data.frame(cell_id=b$cell_ids, truth=truth, cluster=cluster),
          file.path("outputs", dataset, "cluster_labels.csv"), row.names=FALSE)
print(metrics)
