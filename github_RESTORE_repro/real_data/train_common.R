source("R/model_formula.R")

project_simplex <- function(x) {
  u <- sort(x, decreasing=TRUE); cssv <- cumsum(u)-1
  rho <- max(which(u-cssv/seq_along(u)>0)); pmax(x-cssv[rho]/rho,0)
}
normalize_rows <- function(x,target) x/sqrt(rowSums(x^2)+1e-8)*target

train_tmdsp <- function(dataset) {
  cfg <- read.csv("config/datasets.csv", stringsAsFactors=FALSE)
  p <- cfg[cfg$dataset==dataset,,drop=FALSE]
  stem <- if(dataset=="lung") "LG_ACCPU_5kb_573cell" else "M1C_H1930001_5kb_542cell"
  cov <- as.matrix(readRDS(file.path("data/real",paste0(stem,"_cov.rds"))))
  meth <- as.matrix(readRDS(file.path("data/real",paste0(stem,"_methy.rds"))))
  cov[is.na(cov)] <- 0; meth[is.na(meth)] <- 0; mask <- cov>0
  y <- matrix(0,nrow(cov),ncol(cov)); y[mask] <- meth[mask]/cov[mask]
  cell_ids <- colnames(cov); site_ids <- rownames(cov); ns <- nrow(cov); nc <- ncol(cov)

  init_file <- file.path("data/initialization",paste0(dataset,"_state.rds"))
  init <- readRDS(init_file)
  stopifnot(identical(init$cell_ids,cell_ids), identical(init$site_ids,site_ids),
            init$k_rank==p$k_rank)
  U <- init$U_hat; s <- init$s_hat; M <- init$M_hat; R <- init$R_mat
  rownames(U)<-site_ids; rownames(R)<-cell_ids

  W <- Matrix::bandSparse(ns,k=1,diagonals=list(rep(1,ns-1)),symmetric=TRUE)
  L <- Matrix::Diagonal(ns)-Matrix::Diagonal(x=1/(Matrix::rowSums(W)+1e-8))%*%W
  site_cov <- rowMeans(cov); cell_cov <- colMeans(cov)
  cs <- outer(site_cov/(median(site_cov[site_cov>0])+1e-8),
              cell_cov/(median(cell_cov[cell_cov>0])+1e-8),"+")
  coverage_factor <- cs/(cs+1)
  support <- pmin(1,rowSums(mask)/(median(rowSums(mask)[rowSums(mask)>0])+1e-8))
  support <- matrix(support,ns,nc)
  Z <- matrix(.5,ns,nc); dual <- matrix(0,ns,nc)
  kappa0 <- .8; lambda <- p$lambda; main_prob <- .75
  lambda3 <- 1.5; lrU <- .005; lrM <- .005; lrR <- .0005
  max_iter <- as.integer(Sys.getenv("TMDSP_MAX_ITER",p$max_iter))
  for(iter in seq_len(max_iter)) {
    V <- R%*%M; eta <- pmin(pmax(U%*%t(V)+matrix(s,ns,nc),-12),12); prob <- sigmoid(eta)
    adaptive_weight <- pmax(0,1-4*prob*(1-prob))*coverage_factor*support
    kij <- kappa_lambda(kappa0,lambda,adaptive_weight)
    grad <- (prob-y)*mask + kij*(prob-Z+dual)
    U <- U-lrU*pmin(pmax(grad%*%V,-10),10)
    U <- normalize_rows(scale(U,center=TRUE,scale=FALSE),sqrt(p$k_rank))
    gv <- t(grad)%*%U
    M <- M-lrM*pmin(pmax(t(R)%*%gv+.01*M,-10),10)
    M <- normalize_rows(scale(M,center=TRUE,scale=FALSE),sqrt(p$k_rank))
    R <- t(apply(R-lrR*pmin(pmax(gv%*%t(M),-10),10),1,project_simplex))
    gs <- rowSums(grad)+lambda3*as.numeric(L%*%s)+.05*s
    s <- pmin(pmax(s-(lrU*.2)*pmin(pmax(gs,-1),1),-8),8)
    xt <- prob+dual; Z <- if(iter>100) ifelse(xt<.02,0,ifelse(xt>.95,1,xt)) else prob
    Z <- pmin(pmax(Z,0),1); dual <- pmin(pmax(dual+.1*(prob-Z),-.5),.5)
    if(iter%%20==0) message(dataset," iteration ",iter,"/",max_iter)
  }
  V <- R%*%M; fixed <- U%*%t(V)
  shift <- rowSums(fixed*mask)/(rowSums(mask)+1e-8); fixed <- fixed-matrix(shift,ns,nc); s <- s+shift
  for(di in seq_len(30L)) {
    prob <- sigmoid(pmin(pmax(fixed+matrix(s,ns,nc),-12),12))
    step <- rowSums((prob-y)*mask)/(rowSums(mask*prob*(1-prob))+1e-8)
    step <- pmin(pmax(step,-1),1); s <- pmin(pmax(s-step,-8),8)
    if(mean(abs(step))<1e-5) break
  }
  shift <- rowSums(fixed*mask)/(rowSums(mask)+1e-8); fixed <- fixed-matrix(shift,ns,nc); s <- s+shift
  P <- sigmoid(pmin(pmax(fixed+matrix(s,ns,nc),-12),12))
  rownames(P)<-site_ids; colnames(P)<-cell_ids; rownames(V)<-cell_ids
  out <- file.path("outputs",dataset); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(P_hat=P,V_hat=V,cell_ids=cell_ids,lambda=lambda,
    formula="kappa0*((1-lambda)+lambda*w_ij)"),file.path(out,"model.rds"))
  write.csv(data.frame(dataset=dataset,lambda=lambda,kappa0=kappa0,k_rank=p$k_rank,
    max_iter=max_iter,initialization_file=basename(init_file)),
    file.path(out,"training_parameters.csv"),row.names=FALSE)
}
