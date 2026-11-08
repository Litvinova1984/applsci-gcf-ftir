# Resampling-seed reconciliation for Applied Sciences major revision.
# Purpose: compare original R/03_cluster_stability.R results with the
# new variant sensitivity analyses, using the SAME RNG seed for draw pools.
# Does NOT edit the frozen 15 October release or the current D outputs.
# Run from the root: source('R/09b_resampling_seed_audit.R')

# Additional analyses requested during major revision of applsci-4662734.
# Run FROM THE ROOT of the FROZEN 15 October 2026 exploratory release.
# source('R/09_major_revision_sensitivity.R')
# NO baseline files or published tables are overwritten.
# Requires mclust, cluster, FactoMineR (already used by the release).

if (!file.exists(file.path('R', '00_utils.R'))) {
  stop('Run from the root of the October 15 exploratory repository.')
}
source(file.path('R', '00_utils.R'))
for (pkg in c('mclust', 'cluster', 'FactoMineR')) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Missing original dependency: ', pkg)
}

OUT <- file.path('major_revision_outputs','seed_aligned')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
P <- file.path('data', 'pediatric', 'processed')
required <- c(file.path(P, 'pediatric_X_snv.rds'),
              file.path(P, 'pediatric_X_active_absorbance.rds'),
              file.path(P, 'pediatric_wavenumbers.rds'),
              file.path(P, 'pediatric_clusters_baseline.rds'),
              file.path(P, 'pediatric_PCA_scores_3PC.rds'),
              file.path('data','pediatric','pediatric_metadata.csv'),
              file.path('data', 'technical_controls', sprintf('etalon_%02d.dpt', 1:4)))
miss <- required[!file.exists(required)]
if (length(miss)) stop('Missing INPUTS from the frozen release: \n', paste(miss, collapse='\n'))

X <- as.matrix(readRDS(file.path(P, 'pediatric_X_snv.rds')))
Xabs <- as.matrix(readRDS(file.path(P, 'pediatric_X_active_absorbance.rds')))
wn <- as.numeric(readRDS(file.path(P, 'pediatric_wavenumbers.rds')))
meta <- utils::read.csv(file.path('data','pediatric','pediatric_metadata.csv'),
                        stringsAsFactors = FALSE, check.names = FALSE)
base_raw <- readRDS(file.path(P, 'pediatric_clusters_baseline.rds'))
score_frozen <- as.matrix(readRDS(file.path(P, 'pediatric_PCA_scores_3PC.rds')))
stopifnot(identical(dim(X), c(25L,793L)), identical(dim(Xabs), c(25L,793L)),
          length(wn)==793, all(is.finite(X)), all(is.finite(Xabs)))
ids <- rownames(X)
if (is.null(ids) || !setequal(ids, as.character(meta$participant_id))) stop('ID mismatch')
if (is.null(rownames(Xabs)) || !setequal(rownames(Xabs), ids)) stop('Absorbance ID mismatch')
Xabs <- Xabs[ids,,drop=FALSE]
if (!identical(colnames(Xabs), colnames(X))) stop('Absorbance/SNV variable ordering mismatch')
meta <- meta[match(ids, as.character(meta$participant_id)),,drop=FALSE]
meta$sex <- factor(meta$sex, levels=c('F','M'))
meta$age <- as.numeric(meta$age)
if (anyNA(meta$sex) || anyNA(meta$age)) stop('Sex/age missing')
if (is.null(names(base_raw))) {
  base <- as.integer(base_raw)
} else {
  base <- as.integer(as.character(base_raw[ids]))
}
names(base) <- ids
if (!identical(as.integer(table(base)), c(15L,4L,6L))) stop('Unexpected baseline cluster sizes')
if (!is.null(rownames(score_frozen))) score_frozen <- score_frozen[ids,,drop=FALSE]
stopifnot(nrow(score_frozen)==25, ncol(score_frozen)>=3)

write_csv <- function(z, filename) {
  utils::write.csv(z, file=file.path(OUT, filename), row.names=FALSE, na='')
}
center_cols <- function(M) sweep(M, 2, colMeans(M), '-')
center_rows <- function(M) sweep(M, 1, rowMeans(M), '-')
ssq <- function(z) sum(z*z)
fmt <- function(x) sprintf('%.6f', x)

# Low-dimensional Gram PCA yields the same Euclidean PC-score geometry
# as FactoMineR::PCA(X, scale.unit=FALSE). Signs are immaterial for clustering.
fit_partition <- function(M, k=3L, n_pc=3L, return_axes=FALSE) {
  M <- as.matrix(M)
  mu <- colMeans(M)
  C <- sweep(M, 2, mu, '-')
  ev <- eigen(tcrossprod(C), symmetric=TRUE)
  lam <- pmax(ev$values[seq_len(n_pc)], 0)
  if (any(lam <= 1e-10)) stop('Too few nonzero PCs in a resample')
  U <- ev$vectors[,seq_len(n_pc),drop=FALSE]
  scores <- sweep(U, 2, sqrt(lam), '*')
  d <- stats::dist(scores)
  cl <- stats::cutree(stats::hclust(d, method='ward.D2'), k=k)
  out <- list(labels=cl, scores=scores, explained=sum(lam)/sum(pmax(ev$values,0)),
              silhouette=mean(cluster::silhouette(cl,d)[,'sil_width']))
  if (return_axes) {
    V <- crossprod(C,U)
    V <- sweep(V,2,sqrt(lam),'/')
    out$mu <- mu
    out$axes <- V
  }
  out
}
ari <- function(a,b) as.numeric(mclust::adjustedRandIndex(a,b))
sex_residual <- function(M, sex) {
  D <- cbind(intercept=1, sex_m=as.numeric(sex=='M'))
  M - D %*% qr.solve(D,M)
}

# Baseline integrity checks: do not proceed if inputs do not reproduce paper.
primary <- fit_partition(X, k=3L)
base_check <- ari(base, primary$labels)
sex_obs_mat <- sex_residual(X, meta$sex)
sex_obs_ari <- ari(base,fit_partition(sex_obs_mat,k=3L)$labels)
if (base_check < 1 - 1e-10 || abs(sex_obs_ari - 0.134) > 0.003) {
  stop('FROZEN-BASELINE CHECK FAILED: ARI baseline=',fmt(base_check),
       ' and sex residualized=',fmt(sex_obs_ari),
       '. Do not interpret new results until this is resolved.')
}
cat('PASS: frozen baseline 15/4/6; observed sex-residualized ARI=',fmt(sex_obs_ari),'\n')


# Use the canonical seed in R/03_cluster_stability.R:
# set.seed(SEED) before EACH of the subsampling and bootstrap loops.
# This isolates differing random draws as the main likely source of
# small differences between original and initial 09 output.

# D. Alternative specifications.
# Primary fingerprint-only: use ORIGINAL retained absorbance and RECOMPUTE SNV
# within 1800--870 cm^-1. A fixed-original-SNV subset is separately reported.
fidx <- which(wn>=870 & wn<=1800)
if (length(fidx)<200 || length(fidx)>=793) stop('Unexpected fingerprint-only grid')
Xfp <- snv_matrix(Xabs[,fidx,drop=FALSE])
Xfp_sub <- X[,fidx,drop=FALSE]
variants <- list(
  baseline_original_k3=list(M=X,k=3L),
  original_k2=list(M=X,k=2L),
  fingerprint_reSNV_k3=list(M=Xfp,k=3L),
  fingerprint_reSNV_k2=list(M=Xfp,k=2L),
  fingerprint_fixedSNV_k3=list(M=Xfp_sub,k=3L)
)
set.seed(SEED)
sub_indices <- replicate(1000L,sample.int(nrow(X),20L,replace=FALSE),simplify=FALSE)
set.seed(SEED)
boot_indices <- replicate(1000L,sample.int(nrow(X),25L,replace=TRUE),simplify=FALSE)
results <- vector('list',length(variants)); names(results) <- names(variants)

for (nm in names(variants)) {
  spec <- variants[[nm]]
  M <- spec$M; k <- spec$k
  full <- fit_partition(M,k=k)
  labels <- full$labels
  loo <- numeric(nrow(M))
  for(i in seq_len(nrow(M))){
    keep <- setdiff(seq_len(nrow(M)),i)
    loo[i] <- ari(labels[keep],fit_partition(M[keep,,drop=FALSE],k=k)$labels)
  }
  ss <- vapply(sub_indices, function(idx) {
    ari(labels[idx],fit_partition(M[idx,,drop=FALSE],k=k)$labels)
  },numeric(1))
  boot <- vapply(boot_indices, function(idx) {
    fi <- fit_partition(M[idx,,drop=FALSE],k=k,return_axes=TRUE)
    centroids <- do.call(rbind,lapply(sort(unique(fi$labels)),function(g){
      colMeans(fi$scores[fi$labels==g,,drop=FALSE])
    }))
    orig_scores <- sweep(M,2,fi$mu,'-') %*% fi$axes
    ds <- sapply(seq_len(nrow(centroids)),function(j){
      rowSums(sweep(orig_scores,2,centroids[j,],'-')^2)
    })
    if(is.null(dim(ds))) ds <- matrix(ds,nrow=nrow(M))
    projected <- apply(ds,1,which.min)
    ari(labels,projected)
  },numeric(1))
  results[[nm]] <- data.frame(specification=nm,n_wavenumbers=ncol(M),k=k,
      sizes=paste(as.integer(table(labels)),collapse='/'),
      ARI_vs_frozen_k3=ari(base,labels),
      PC123_variance_percent=100*full$explained,
      silhouette=full$silhouette,
      LOO_mean_ARI=mean(loo),LOO_median_ARI=median(loo),
      subsample_mean_ARI=mean(ss),subsample_median_ARI=median(ss),
      bootstrap_mean_ARI=mean(boot),bootstrap_median_ARI=median(boot))
  cat('D complete:',nm,' | ARI vs original=',fmt(ari(base,labels)),'\n')
}
write_csv(do.call(rbind,results),'D_alternative_windows_and_k.csv')


# Check differences from archived October 15 summary (do not overwrite it).
all_out <- do.call(rbind, results)
row <- all_out[all_out$specification=='baseline_original_k3',,drop=FALSE]
original_subsample <- 0.6418108439
original_bootstrap <- 0.5794283128
cat('\n=== Resampling alignment audit ===\n')
cat('Original subsampling mean ARI: ', sprintf('%.10f',original_subsample),'\n',sep='')
cat('Aligned subsampling mean ARI:  ',sprintf('%.10f',row$subsample_mean_ARI),'\n',sep='')
cat('Original bootstrap mean ARI:   ',sprintf('%.10f',original_bootstrap),'\n',sep='')
cat('Aligned bootstrap mean ARI:    ',sprintf('%.10f',row$bootstrap_mean_ARI),'\n',sep='')

comparison <- data.frame(
  metric=c('subsample_mean_ARI','bootstrap_mean_ARI'),
  manuscript=c(original_subsample, original_bootstrap),
  seed_aligned=c(row$subsample_mean_ARI,row$bootstrap_mean_ARI)
)
comparison$absolute_difference <- abs(comparison$seed_aligned - comparison$manuscript)
write_csv(comparison,'09b_resampling_vs_manuscript.csv')
if (all(comparison$absolute_difference < 0.000001)) {
  cat('PASS: both means match the original report to 1e-6.\n')
} else {
  cat('REVIEW: differences remain. Check sample ordering, PCA/centroid algorithm, and R RNG kind.\n')
}
cat('Outputs saved in: ', normalizePath(OUT),'\n',sep='')
