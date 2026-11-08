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

OUT <- file.path('major_revision_outputs')
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

# A. Requested sex-label permutation benchmark: group sizes 12 F / 13 M.
B <- 4999L
set.seed(SEED + 901L)
perm_ari <- numeric(B)
for (b in seq_len(B)) {
  perm_sex <- sample(meta$sex, replace=FALSE)
  perm_ari[b] <- ari(base,fit_partition(sex_residual(X,perm_sex),k=3L)$labels)
  if (b %% 1000L == 0L) cat('Permutation',b,'/',B,'\n')
}
perm_p <- (1 + sum(perm_ari <= sex_obs_ari + 1e-12)) / (B+1)
write_csv(data.frame(permutation=seq_len(B),ARI_vs_baseline=perm_ari),
          'A_sex_permutations_4999.csv')
write_csv(data.frame(observed_ARI=sex_obs_ari, B=B,
                     lower_tail_permutation_p=perm_p,
                     permutation_median=median(perm_ari),
                     permutation_q025=unname(quantile(perm_ari,.025)),
                     permutation_q975=unname(quantile(perm_ari,.975))),
          'A_sex_permutation_summary.csv')

# B. Variance removed by sex or clean-strip direction.
# Reconstruct the reference EXACTLY as in R/05_strip_reproducibility.R:
# individual strip -> SNV -> mean of the four -> center across wavenumbers -> unit vector.
read_strip <- function(f) {
  dat <- utils::read.csv(f,header=FALSE,stringsAsFactors=FALSE)
  if (ncol(dat)!=2L) stop('Expected 2-column DPT: ',f)
  dat <- dat[is.finite(dat[[1]]) & is.finite(dat[[2]]),,drop=FALSE]
  dat <- dat[(dat[[1]]>=2800 & dat[[1]]<=3400) |
             (dat[[1]]>=870 & dat[[1]]<=1800),,drop=FALSE]
  dat <- dat[order(dat[[1]],decreasing=TRUE),,drop=FALSE]
  if (nrow(dat)!=length(wn) || max(abs(dat[[1]]-wn))>1e-4) {
    stop('Technical-strip wavenumber grid differs from frozen pediatric grid: ',f)
  }
  as.numeric(dat[[2]])
}
strip_files <- file.path('data','technical_controls',sprintf('etalon_%02d.dpt',1:4))
strip_raw <- do.call(rbind,lapply(strip_files,read_strip))
strip_snv <- t(apply(strip_raw,1,snv))
reference <- colMeans(strip_snv)
e <- reference - mean(reference)
e <- e/sqrt(ssq(e))
Xrowc <- center_rows(X)
coeff <- as.numeric(Xrowc %*% e)
Xclean <- Xrowc - tcrossprod(coeff,e)
clean_ari <- ari(base,fit_partition(Xclean,k=3L)$labels)
if (abs(clean_ari-0.407)>0.003) {
  stop('CLEAN-STRIP CHECK FAILED: ARI=',fmt(clean_ari),
       '; frozen value 0.407. Check strip alignment/preprocessing.')
}
Xc <- center_cols(X)
XsexC <- center_cols(sex_obs_mat)
XcleanC <- center_cols(Xclean)
Rsex <- Xc-XsexC
Rclean <- Xc-XcleanC
frac_sex <- ssq(Rsex)/ssq(Xc)
frac_clean <- ssq(Rclean)/ssq(Xc)
cosine_removed <- sum(Rsex*Rclean)/sqrt(ssq(Rsex)*ssq(Rclean))
# Relative changes to between-PARTICIPANT spectral variance after column centering.
# This is NOT a quantitative biological/substrate partition.
write_csv(data.frame(sex_removed_percent=100*frac_sex,
                     clean_direction_removed_percent=100*frac_clean,
                     frobenius_cosine=cosine_removed,
                     squared_cosine=cosine_removed^2,
                     clean_coefficient_by_sex_cor=cor(coeff,as.numeric(meta$sex=='M')),
                     sex_ARI=sex_obs_ari,clean_ARI=clean_ari),
          'B_removed_variation_and_overlap.csv')
write_csv(data.frame(participant_id=ids,sex=as.character(meta$sex),
                     clean_projection_coefficient=coeff,
                     clean_component_energy=coeff^2),
          'B_clean_projection_by_participant.csv')

# C. Age and available puberty by sex. The existing PCA scores are used
# to preserve the frozen sign convention of PC1.
sex_level <- c('F','M')
age_by <- do.call(rbind,lapply(sex_level,function(s){
  x <- meta$age[meta$sex==s]
  data.frame(sex=s,n=length(x),age_mean=mean(x),age_SD=sd(x),
             age_median=median(x),age_min=min(x),age_max=max(x))
}))
write_csv(age_by,'C_age_by_sex.csv')
pub_levels <- c('Normal','Early','Delayed','Unknown')
pub <- as.data.frame(table(sex=factor(meta$sex,levels=sex_level),
                           puberty=factor(meta$pubertal_development,levels=pub_levels)))
write_csv(pub,'C_puberty_by_sex.csv')
meta$PC1 <- as.numeric(score_frozen[,1])
fit0 <- lm(PC1 ~ sex,data=meta)
fit1 <- lm(PC1 ~ sex + age,data=meta)
lm_tab <- function(fit,label){
  z <- as.data.frame(summary(fit)$coefficients)
  ci <- confint(fit)
  data.frame(model=label,term=rownames(z),estimate=z[[1]],SE=z[[2]],
             t=z[[3]],p_value=z[[4]],ci_low=ci[,1],ci_high=ci[,2],
             row.names=NULL,check.names=FALSE)
}
write_csv(rbind(lm_tab(fit0,'PC1 ~ sex'),lm_tab(fit1,'PC1 ~ sex + age')),
          'C_PC1_sex_age_linear_sensitivity.csv')
# HC3 (heteroscedasticity-robust) standard errors without additional packages.
D_hc3 <- model.matrix(fit1)
h <- rowSums((D_hc3 %*% solve(crossprod(D_hc3))) * D_hc3)
u_hc3 <- resid(fit1)/(1-h)
V_hc3 <- solve(crossprod(D_hc3)) %*% crossprod(D_hc3 * as.numeric(u_hc3)) %*%
         solve(crossprod(D_hc3))
se_hc3 <- sqrt(diag(V_hc3))
t_hc3 <- coef(fit1)/se_hc3
df_res <- df.residual(fit1)
write_csv(data.frame(term=names(coef(fit1)),estimate=as.numeric(coef(fit1)),
                     HC3_SE=se_hc3,HC3_t=t_hc3,
                     HC3_p=2*pt(-abs(t_hc3),df=df_res),
                     HC3_CI_low=coef(fit1)-qt(.975,df_res)*se_hc3,
                     HC3_CI_high=coef(fit1)+qt(.975,df_res)*se_hc3),
          'C_PC1_sex_age_HC3_robust.csv')
write_csv(as.data.frame(table(sex=meta$sex,endocrine=meta$endocrine_status)),
          'C_endocrine_status_by_sex_descriptive.csv')
write_csv(as.data.frame(table(sex=meta$sex,uctd=meta$uctd_status)),
          'C_UCTD_status_by_sex_descriptive.csv')

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
set.seed(SEED + 1900L)
sub_indices <- replicate(1000L,sample.int(nrow(X),20L,replace=FALSE),simplify=FALSE)
set.seed(SEED + 1901L)
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

# Output manifest and a compact review summary.
files <- list.files(OUT,pattern='\\.csv$',full.names=FALSE)
write_csv(data.frame(output=files), '00_major_revision_output_manifest.csv')
cat('\nSUCCESS: outputs saved under ',normalizePath(OUT),'\n',sep='')
cat('Key values: sex-permutation lower-tail p=',fmt(perm_p),
    '; removed variation sex=',fmt(100*frac_sex),'%; clean=',
    fmt(100*frac_clean),'%.\n')
cat('Baseline integrity checked: sex ARI 0.134; clean ARI 0.407.\n')
