# Applied Sciences applsci-4662734 -- Major Revision
# Reproduce the PRELIMINARY analyses already reported in the response:
#   (C) age and available puberty status by sex;
#   (C) PC1 ~ sex and PC1 ~ sex + age (standard and HC3 inference);
#   (D) k=2 vs k=3 using the *frozen* PC1--PC3 scores.
#
# Reproducibility basis: v1.0.0-pediatric-exploratory (15 October 2026)
# https://github.com/Litvinova1984/applsci-gcf-ftir/releases/tag/v1.0.0-pediatric-exploratory
#
# IMPORTANT:
# - Run from the repository ROOT in RStudio:
#     source("R/09a_major_revision_preliminary_reproduction.R")
# - No refitting of PCA is performed here. No original published result is overwritten.
# - The k=2 calculation uses the published/frozen 3-PC scores and is NOT a
#   resampling stability test; that belongs to 09_major_revision_sensitivity.R.
# - No new biological or causal inference is implied by the age adjustment.
# - This file was prepared without an R runtime and MUST be run and checked
#   against the frozen release before uploading numerical results to GitHub.

DATA_DIR <- file.path("data", "pediatric")
OBJ_DIR <- file.path(DATA_DIR, "processed")
OUTPUT_DIR <- file.path("major_revision_outputs", "preliminary")
required_files <- c(
  metadata = file.path(DATA_DIR, "pediatric_metadata.csv"),
  frozen_PCA_scores = file.path(OBJ_DIR, "pediatric_PCA_scores_3PC.rds"),
  frozen_k3_labels = file.path(OBJ_DIR, "pediatric_clusters_baseline.rds")
)
not_found <- required_files[!file.exists(required_files)]
if (length(not_found)) {
  stop("Required frozen-release inputs missing:\n  ",
       paste(not_found, collapse = "\n  "),
       "\nRun from the 15 October 2026 exploratory repository root.")
}
if (!requireNamespace("cluster", quietly = TRUE)) {
  stop("The 'cluster' package is required (same silhouette implementation as the article).")
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_output <- function(dat, name) {
  utils::write.csv(dat, file.path(OUTPUT_DIR, name), row.names = FALSE,
                   na = "", fileEncoding = "UTF-8")
}
near <- function(x, expected, tol = 1e-4, label = "quantity") {
  if (length(x) != 1L || !is.finite(x) || abs(x - expected) > tol) {
    stop(label, " did not reproduce the previously reported number. Obtained: ",
         format(x, digits = 12), "; expected: ",
         format(expected, digits = 12), ". STOP and check input version/order.")
  }
}

# I. Read and explicitly align frozen inputs -----------------------------------
meta <- utils::read.csv(required_files[["metadata"]],
                        stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("participant_id", "sex", "age", "pubertal_development")
if (!all(needed %in% colnames(meta))) {
  stop("pediatric_metadata.csv is missing: ",
       paste(setdiff(needed, names(meta)), collapse = ", "))
}
meta$participant_id <- as.character(meta$participant_id)
if (nrow(meta) != 25L || anyDuplicated(meta$participant_id)) {
  stop("Expected exactly 25 distinct participant IDs.")
}
meta$sex <- factor(meta$sex, levels = c("F", "M"))
meta$age <- as.numeric(meta$age)
if (anyNA(meta$sex) || anyNA(meta$age)) stop("Missing or invalid age/sex.")
if (!identical(as.integer(table(meta$sex)), c(12L, 13L))) {
  stop("Unexpected sex counts; this does not match the frozen cohort.")
}

score_object <- readRDS(required_files[["frozen_PCA_scores"]])
if (!is.matrix(score_object) && !is.data.frame(score_object)) {
  stop("Unexpected PCA scores object type; inspect the frozen file.")
}
score_matrix <- as.matrix(score_object)
if (nrow(score_matrix) != 25L) stop("Expected 25 rows in frozen PCA scores.")
if (all(c("Dim.1", "Dim.2", "Dim.3") %in% colnames(score_matrix))) {
  score_matrix <- score_matrix[, c("Dim.1", "Dim.2", "Dim.3"), drop = FALSE]
} else if (all(c("PC1", "PC2", "PC3") %in% colnames(score_matrix))) {
  score_matrix <- score_matrix[, c("PC1", "PC2", "PC3"), drop = FALSE]
} else if (ncol(score_matrix) == 3L) {
  # The released object consists only of three PC-score columns.
  score_matrix <- score_matrix[, 1:3, drop = FALSE]
} else {
  stop("Cannot identify the three PCA-score columns in the frozen object.")
}
storage.mode(score_matrix) <- "double"
if (any(!is.finite(score_matrix))) stop("Nonfinite frozen PC scores.")
colnames(score_matrix) <- c("PC1", "PC2", "PC3")

# If score rows have IDs, align clinical metadata by them. If not, require the
# explicit 1:25 order in the frozen metadata instead of silently guessing.
score_ids <- rownames(score_matrix)
if (!is.null(score_ids) && setequal(score_ids, meta$participant_id)) {
  meta <- meta[match(score_ids, meta$participant_id), , drop = FALSE]
} else if (is.null(score_ids) || identical(score_ids, as.character(1:25))) {
  if (!identical(meta$participant_id, as.character(1:25))) {
    stop("Unlabelled PCA rows; metadata order is not 1:25 -- unsafe alignment.")
  }
} else {
  stop("Frozen PCA row names do not match participant IDs. Inspect manually.")
}
stopifnot(nrow(meta) == nrow(score_matrix))

raw_labels <- readRDS(required_files[["frozen_k3_labels"]])
if (is.matrix(raw_labels) && ncol(raw_labels) == 1L) raw_labels <- as.vector(raw_labels[, 1L])
if (!is.atomic(raw_labels) || length(raw_labels) != 25L) {
  stop("Unexpected frozen baseline label format.")
}
if (!is.null(names(raw_labels)) && setequal(names(raw_labels), meta$participant_id)) {
  raw_labels <- raw_labels[match(meta$participant_id, names(raw_labels))]
} else if (!is.null(names(raw_labels)) &&
           !identical(names(raw_labels), as.character(1:25))) {
  stop("Frozen baseline labels have unmatched ID names; unsafe to proceed.")
}
k3_frozen <- as.character(raw_labels)
if (!identical(sort(as.integer(table(k3_frozen))), c(4L, 6L, 15L))) {
  stop("Baseline cluster sizes do not match 15 / 4 / 6.")
}

# II. Descriptive age and recorded pubertal category by sex --------------------
sex_groups <- c("F", "M")
age_table <- do.call(rbind, lapply(sex_groups, function(g) {
  x <- meta$age[meta$sex == g]
  data.frame(sex = g, n = length(x), mean_age = mean(x), SD_age = sd(x),
             median_age = median(x), min_age = min(x), max_age = max(x))
}))
write_output(age_table, "C_age_by_sex.csv")

puberty_levels <- c("Normal", "Early", "Delayed", "Unknown")
if (!all(meta$pubertal_development %in% puberty_levels)) {
  stop("Unexpected puberty category; preserve the original labels explicitly.")
}
puberty_tab <- table(sex = factor(meta$sex, levels = sex_groups),
                     pubertal_development = factor(meta$pubertal_development,
                                                   levels = puberty_levels))
puberty_wide <- data.frame(sex = rownames(puberty_tab),
                           as.data.frame.matrix(puberty_tab), row.names = NULL,
                           check.names = FALSE)
write_output(puberty_wide, "C_puberty_by_sex.csv")
near(age_table$mean_age[age_table$sex == "F"], 11.4166666667, label = "Mean F age")
near(age_table$mean_age[age_table$sex == "M"], 12.3846153846, label = "Mean M age")
if (!identical(as.integer(puberty_tab["F", ]), c(4L, 3L, 0L, 5L)) ||
    !identical(as.integer(puberty_tab["M", ]), c(4L, 0L, 1L, 8L))) {
  stop("Pubertal category counts differ from the previously reported table.")
}

# III. PC1 sex association and age-adjusted descriptive sensitivity ------------
meta$PC1 <- as.numeric(score_matrix[, 1L])
fit_sex <- stats::lm(PC1 ~ sex, data = meta)
fit_age <- stats::lm(PC1 ~ sex + age, data = meta)
model_rows <- function(fit, model_name) {
  tab <- summary(fit)$coefficients
  ci <- stats::confint(fit, level = 0.95)
  data.frame(model = model_name, term = rownames(tab),
             estimate = unname(tab[, "Estimate"]),
             standard_error = unname(tab[, "Std. Error"]),
             t_value = unname(tab[, "t value"]),
             p_value = unname(tab[, "Pr(>|t|)"]),
             CI_low = unname(ci[, 1L]), CI_high = unname(ci[, 2L]),
             row.names = NULL, check.names = FALSE)
}
write_output(rbind(model_rows(fit_sex, "PC1 ~ sex"),
                   model_rows(fit_age, "PC1 ~ sex + age")),
             "C_PC1_sex_age_standard.csv")

# HC3 heteroscedasticity-robust covariance, with t(df=n-p) inference.
Z <- stats::model.matrix(fit_age)
Z_inv <- solve(crossprod(Z))
leverages <- rowSums((Z %*% Z_inv) * Z)
if (any(leverages >= 1)) stop("HC3: leverage 1; cannot compute covariance.")
hc3_weighted <- sweep(Z, 1L, stats::residuals(fit_age) / (1 - leverages), "*")
hc3_cov <- Z_inv %*% crossprod(hc3_weighted) %*% Z_inv
hc3_SE <- sqrt(diag(hc3_cov))
hc3_t <- as.numeric(stats::coef(fit_age)) / hc3_SE
hc3_df <- stats::df.residual(fit_age)
hc3_crit <- stats::qt(0.975, df = hc3_df)
hc3_table <- data.frame(term = names(stats::coef(fit_age)),
                        estimate = unname(stats::coef(fit_age)),
                        HC3_SE = unname(hc3_SE),
                        HC3_t = hc3_t,
                        HC3_p = 2 * stats::pt(-abs(hc3_t), df = hc3_df),
                        HC3_CI_low = unname(stats::coef(fit_age) - hc3_crit * hc3_SE),
                        HC3_CI_high = unname(stats::coef(fit_age) + hc3_crit * hc3_SE),
                        row.names = NULL)
write_output(hc3_table, "C_PC1_sex_age_HC3.csv")

sex_est_unadj <- unname(stats::coef(fit_sex)["sexM"])
sex_est_adj <- unname(stats::coef(fit_age)["sexM"])
sex_p_adj <- unname(summary(fit_age)$coefficients["sexM", "Pr(>|t|)"])
sex_hc3_p <- unname(hc3_table$HC3_p[hc3_table$term == "sexM"])
near(sex_est_unadj, 6.8747211202, 0.002, "PC1 male-minus-female, unadjusted")
near(sex_est_adj, 6.7539409262, 0.002, "PC1 male-minus-female, age-adjusted")
near(sex_p_adj, 0.0084737571, 0.0001, "PC1 age-adjusted classical p")
near(sex_hc3_p, 0.0161, 0.001, "PC1 age-adjusted HC3 p")

# IV. Recalculate k=2 and k=3 on the SAME frozen PC-score geometry -------------
# Only comparisons of frozen scores; not a new PCA fit, not LOAO / resampling.
D <- stats::dist(score_matrix, method = "euclidean")
H <- stats::hclust(D, method = "ward.D2")
k2 <- stats::cutree(H, k = 2L)
k3_new <- stats::cutree(H, k = 3L)

# Adjusted Rand Index -- implemented in base R, invariant to label names.
choose2 <- function(v) sum(as.numeric(v) * (as.numeric(v) - 1) / 2)
adjusted_rand <- function(a, b) {
  tab <- table(as.character(a), as.character(b))
  n <- sum(tab)
  observed <- choose2(tab)
  marginal_row <- choose2(rowSums(tab))
  marginal_col <- choose2(colSums(tab))
  expected <- marginal_row * marginal_col / choose2(n)
  max_index <- 0.5 * (marginal_row + marginal_col)
  if (abs(max_index - expected) < 1e-12) {
    return(if (all(a == b)) 1 else 0)
  }
  (observed - expected) / (max_index - expected)
}

# Calinski--Harabasz (between / within cluster sums of squares).
calinski_harabasz <- function(coords, labels) {
  labs <- unique(as.character(labels))
  n <- nrow(coords); k <- length(labs)
  center <- colMeans(coords)
  total_ss <- sum(sweep(coords, 2L, center, "-")^2)
  within <- sum(vapply(labs, function(g) {
    M <- coords[labels == g, , drop = FALSE]
    sum(sweep(M, 2L, colMeans(M), "-")^2)
  }, numeric(1L)))
  ((total_ss - within) / (k - 1)) / (within / (n - k))
}
# Davies--Bouldin (mean within-cluster Euclidean radius / centroid separation).
davies_bouldin <- function(coords, labels) {
  labs <- unique(as.character(labels)); k <- length(labs)
  centroids <- do.call(rbind, lapply(labs, function(g) {
    colMeans(coords[labels == g, , drop = FALSE])
  }))
  scatter <- vapply(seq_along(labs), function(j) {
    M <- coords[labels == labs[j], , drop = FALSE]
    sqrt_d <- sqrt(rowSums(sweep(M, 2L, centroids[j, ], "-")^2))
    mean(sqrt_d)
  }, numeric(1L))
  r <- matrix(-Inf, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) if (i != j) {
    dij <- sqrt(sum((centroids[i, ] - centroids[j, ])^2))
    r[i, j] <- (scatter[i] + scatter[j]) / dij
  }
  mean(apply(r, 1L, max))
}
cluster_metrics <- function(labels, k) {
  data.frame(k = k,
             group_sizes_sorted = paste(sort(as.integer(table(labels))), collapse = "/"),
             silhouette = mean(cluster::silhouette(labels, D)[, "sil_width"]),
             Calinski_Harabasz = calinski_harabasz(score_matrix, labels),
             Davies_Bouldin = davies_bouldin(score_matrix, labels),
             ARI_vs_frozen_k3 = adjusted_rand(labels, k3_frozen))
}
result_k2 <- cluster_metrics(k2, 2L)
result_k3 <- cluster_metrics(k3_new, 3L)
write_output(rbind(result_k2, result_k3), "D_k2_vs_k3_frozen_PC_scores.csv")
write_output(data.frame(participant_id = meta$participant_id,
                        frozen_k3 = k3_frozen,
                        recalculated_k3 = k3_new,
                        k2 = k2),
             "D_k2_vs_k3_participant_assignments.csv")
near(result_k3$ARI_vs_frozen_k3, 1, 1e-8, "k=3 labels vs frozen labels ARI")
near(result_k3$silhouette, 0.4342488043, 1e-4, "k=3 silhouette")
near(result_k2$ARI_vs_frozen_k3, 0.6147919877, 1e-4, "k=2 vs frozen k=3 ARI")
near(result_k2$silhouette, 0.4376790349, 1e-4, "k=2 silhouette")
near(result_k2$Calinski_Harabasz, 17.6463041087, 1e-3, "k=2 Calinski-Harabasz")
near(result_k2$Davies_Bouldin, 0.8975678104, 1e-4, "k=2 Davies-Bouldin")

# V. Reproducibility manifest and session --------------------------------------
manifest <- data.frame(input = names(required_files),
                       path = unname(required_files),
                       MD5 = unname(tools::md5sum(required_files)))
write_output(manifest, "00_input_MD5_manifest.csv")
writeLines(capture.output(utils::sessionInfo()),
           file.path(OUTPUT_DIR, "00_sessionInfo.txt"))
cat("\nPASS -- preliminary checks reproduced and saved to:",
    normalizePath(OUTPUT_DIR), "\n")
cat("  Age-adjusted M-F PC1:", signif(sex_est_adj, 6),
    ", p:", signif(sex_p_adj, 6),
    ", HC3 p:", signif(sex_hc3_p, 6), "\n")
cat("  k=2 vs frozen k=3 ARI:", signif(result_k2$ARI_vs_frozen_k3, 6),
    ", silhouette:", signif(result_k2$silhouette, 6), "\n")
