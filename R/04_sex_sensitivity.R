# ============================================================
# 04_sex_sensitivity.R
#
# Pediatric GCF FTIR exploratory reframe
# Sex-associated spectral variation and sensitivity analysis
#
# PURPOSE
# -------
# This script retains only the participant-level analysis that is
# central to the reframed Applied Sciences manuscript:
#
#   1. descriptive cohort metadata for Table 1;
#   2. association between sex and the fixed baseline partition;
#   3. sex differences in PC1-PC3, with Benjamini-Hochberg
#      correction across the three planned PC comparisons;
#   4. feature-wise removal of the linear sex-associated component,
#      followed by a fresh PCA and the same Ward.D2 k = 3 clustering;
#   5. comparison of the sex-residualized and baseline partitions.
#
# Age, endocrine status, UCTD status, periodontal status, and
# pubertal-development categories remain in the cohort-description
# table. They are not tested against the data-derived clusters in
# this reframed analysis. This avoids an unfocused set of sparse
# post-clustering clinical comparisons and keeps the inferential
# emphasis on the prespecified sex-sensitivity question.
#
# INTERPRETATION
# --------------
# The residualization is a sensitivity analysis, not a causal model
# or a corrected ground-truth partition. A changed partition shows
# that variation linearly associated with sex contributes to the
# observed multivariate organization in this sample. It does not
# establish a hormonal, developmental, or other biological mechanism.
# ============================================================


# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "00_utils.R"))

required_pkgs <- c(
  "FactoMineR",
  "cluster",
  "mclust",
  "tibble",
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}


# ------------------------------------------------------------
# 1. Paths and fixed analytical settings
# ------------------------------------------------------------

processed_dir <- file.path("data", "pediatric", "processed")
metadata_file <- file.path("data", "pediatric", "pediatric_metadata.csv")
x_file <- file.path(processed_dir, "pediatric_X_snv.rds")
pca_file <- file.path(processed_dir, "pediatric_PCA_baseline_3PC.rds")
scores_file <- file.path(processed_dir, "pediatric_PCA_scores_3PC.rds")
baseline_cluster_file <- file.path(
  processed_dir,
  "pediatric_clusters_baseline.rds"
)

tables_dir <- "tables"
figures_dir <- "figures"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  metadata_file,
  x_file,
  pca_file,
  scores_file,
  baseline_cluster_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) are missing. Run scripts 01 and 02 first:\n",
    paste(missing_files, collapse = "\n")
  )
}

N_PC_BASELINE <- 3L
K_BASELINE <- 3L
N_BOOT_EFFECT <- 5000L


# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

coerce_named_clusters <- function(x, ids, label) {
  if (length(x) != length(ids)) {
    stop(label, " length does not match the participant count.")
  }

  if (is.null(names(x))) {
    out <- if (is.factor(x)) {
      as.integer(as.character(x))
    } else {
      as.integer(x)
    }
    names(out) <- ids
    return(out)
  }

  if (!setequal(names(x), ids)) {
    stop(label, " participant IDs do not match the spectral matrix.")
  }

  x <- x[ids]
  out <- if (is.factor(x)) {
    as.integer(as.character(x))
  } else {
    as.integer(x)
  }
  names(out) <- ids
  out
}

cramers_v <- function(x, y) {
  tab <- table(x, y)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]

  if (nrow(tab) < 2L || ncol(tab) < 2L) {
    return(NA_real_)
  }

  n <- sum(tab)
  expected <- outer(rowSums(tab), colSums(tab)) / n
  valid <- expected > 0
  chi_sq <- sum((tab[valid] - expected[valid])^2 / expected[valid])
  denom <- min(nrow(tab) - 1L, ncol(tab) - 1L)

  sqrt(chi_sq / (n * denom))
}

bootstrap_cramers_v <- function(
    x,
    y,
    B = N_BOOT_EFFECT,
    seed = SEED
) {
  n <- length(x)
  set.seed(seed)

  vals <- replicate(
    B,
    {
      idx <- sample.int(n, size = n, replace = TRUE)
      cramers_v(x[idx], y[idx])
    }
  )

  vals <- vals[is.finite(vals)]

  if (length(vals) == 0L) {
    return(c(lower = NA_real_, upper = NA_real_, valid_bootstrap = 0L))
  }

  c(
    lower = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
    upper = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
    valid_bootstrap = length(vals)
  )
}

rank_biserial_two_group <- function(x_f, x_m) {
  x_f <- x_f[is.finite(x_f)]
  x_m <- x_m[is.finite(x_m)]
  n_f <- length(x_f)
  n_m <- length(x_m)

  if (n_f == 0L || n_m == 0L) {
    return(NA_real_)
  }

  pooled <- c(x_f, x_m)
  ranks <- rank(pooled, ties.method = "average")
  R_f <- sum(ranks[seq_len(n_f)])
  U_f <- R_f - n_f * (n_f + 1) / 2

  2 * U_f / (n_f * n_m) - 1
}

hedges_g_two_group <- function(x_f, x_m) {
  x_f <- x_f[is.finite(x_f)]
  x_m <- x_m[is.finite(x_m)]
  n_f <- length(x_f)
  n_m <- length(x_m)
  df <- n_f + n_m - 2L

  if (n_f < 2L || n_m < 2L || df <= 0L) {
    return(NA_real_)
  }

  pooled_sd <- sqrt(
    ((n_f - 1) * stats::var(x_f) + (n_m - 1) * stats::var(x_m)) / df
  )

  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }

  d <- (mean(x_f) - mean(x_m)) / pooled_sd
  J <- 1 - 3 / (4 * (n_f + n_m) - 9)
  J * d
}

bootstrap_two_group_effect <- function(
    x,
    sex,
    effect_fun,
    B = N_BOOT_EFFECT,
    seed = SEED
) {
  x_f <- x[sex == "F"]
  x_m <- x[sex == "M"]
  set.seed(seed)

  vals <- replicate(
    B,
    {
      boot_f <- sample(x_f, size = length(x_f), replace = TRUE)
      boot_m <- sample(x_m, size = length(x_m), replace = TRUE)
      effect_fun(boot_f, boot_m)
    }
  )

  vals <- vals[is.finite(vals)]

  if (length(vals) == 0L) {
    return(c(lower = NA_real_, upper = NA_real_, valid_bootstrap = 0L))
  }

  c(
    lower = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
    upper = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
    valid_bootstrap = length(vals)
  )
}

ward_clusters <- function(scores, k) {
  hc <- stats::hclust(
    stats::dist(scores, method = "euclidean"),
    method = "ward.D2"
  )
  stats::cutree(hc, k = k)
}

all_permutations <- function(x) {
  if (length(x) == 1L) {
    return(list(x))
  }

  out <- list()
  counter <- 1L

  for (i in seq_along(x)) {
    sub_perm <- all_permutations(x[-i])
    for (p in sub_perm) {
      out[[counter]] <- c(x[i], p)
      counter <- counter + 1L
    }
  }

  out
}

align_cluster_labels <- function(reference, candidate) {
  ref_levels <- sort(unique(reference))
  cand_levels <- sort(unique(candidate))

  if (length(ref_levels) != length(cand_levels)) {
    stop("Reference and candidate partitions have different numbers of clusters.")
  }

  permutations <- all_permutations(ref_levels)
  agreements <- vapply(
    permutations,
    function(p) {
      map <- stats::setNames(p, cand_levels)
      aligned <- unname(map[as.character(candidate)])
      sum(aligned == reference)
    },
    numeric(1)
  )

  best <- which.max(agreements)
  map <- stats::setNames(permutations[[best]], cand_levels)
  aligned <- as.integer(unname(map[as.character(candidate)]))
  names(aligned) <- names(candidate)
  aligned
}

categorical_overall_summary <- function(data, variable) {
  data |>
    dplyr::count(.data[[variable]], name = "n", .drop = FALSE) |>
    dplyr::transmute(
      variable = variable,
      category = as.character(.data[[variable]]),
      n = n,
      mean = NA_real_,
      sd = NA_real_,
      median = NA_real_,
      min = NA_real_,
      max = NA_real_
    )
}


# ------------------------------------------------------------
# 3. Load and validate the frozen analysis cohort
# ------------------------------------------------------------

metadata <- readr::read_csv(metadata_file, show_col_types = FALSE) |>
  dplyr::mutate(
    participant_id = as.character(participant_id),
    age = as.numeric(age),
    sex = factor(sex, levels = c("F", "M")),
    pubertal_development = factor(
      pubertal_development,
      levels = c("Normal", "Early", "Delayed", "Unknown")
    ),
    endocrine_status = factor(
      endocrine_status,
      levels = c(
        "No documented pathology",
        "Documented pathology",
        "Unknown"
      )
    ),
    uctd_status = factor(
      uctd_status,
      levels = c("No documented pathology", "Documented pathology")
    ),
    periodontal_status = factor(
      periodontal_status,
      levels = c("No documented pathology", "Documented pathology")
    )
  )

X_snv <- as.matrix(readRDS(x_file))
pca_baseline <- readRDS(pca_file)
scores_baseline <- as.matrix(readRDS(scores_file))
cluster_baseline_raw <- readRDS(baseline_cluster_file)

sample_ids <- rownames(X_snv)

if (is.null(sample_ids)) {
  stop("X_snv must have participant IDs as row names.")
}

if (!identical(dim(X_snv), c(25L, 793L))) {
  stop(
    "Expected X_snv dimensions 25 x 793; obtained ",
    paste(dim(X_snv), collapse = " x "),
    "."
  )
}

if (any(!is.finite(X_snv))) {
  stop("Non-finite values detected in X_snv.")
}

if (nrow(metadata) != 25L || anyDuplicated(metadata$participant_id) > 0L) {
  stop("The metadata file must contain 25 unique participant IDs.")
}

if (!setequal(metadata$participant_id, sample_ids)) {
  stop("Metadata participant IDs do not match X_snv row names.")
}

metadata <- metadata[
  match(sample_ids, metadata$participant_id),
  ,
  drop = FALSE
]

if (!identical(dim(scores_baseline), c(25L, 3L))) {
  stop("Expected the saved baseline score matrix to be 25 x 3.")
}

if (is.null(rownames(scores_baseline))) {
  rownames(scores_baseline) <- sample_ids
} else if (!setequal(rownames(scores_baseline), sample_ids)) {
  stop("Saved baseline PCA-score IDs do not match X_snv.")
} else {
  scores_baseline <- scores_baseline[sample_ids, , drop = FALSE]
}

cluster_baseline <- coerce_named_clusters(
  cluster_baseline_raw,
  sample_ids,
  "Baseline cluster vector"
)

if (!identical(as.integer(table(cluster_baseline)), c(15L, 4L, 6L))) {
  stop("Expected baseline cluster sizes 15 / 4 / 6.")
}

metadata$baseline_cluster <- factor(
  cluster_baseline[metadata$participant_id],
  levels = 1:K_BASELINE
)

metadata_validation <- tibble::tibble(
  check = c(
    "n_participants",
    "age_min",
    "age_max",
    "age_mean",
    "age_sd",
    "sex_F",
    "sex_M"
  ),
  observed = c(
    nrow(metadata),
    min(metadata$age),
    max(metadata$age),
    mean(metadata$age),
    stats::sd(metadata$age),
    sum(metadata$sex == "F"),
    sum(metadata$sex == "M")
  ),
  expected = c(25, 6, 17, 11.92, 3.37, 12, 13)
) |>
  dplyr::mutate(
    matches = dplyr::case_when(
      check %in% c("age_mean", "age_sd") ~
        abs(round(observed, 2) - expected) < 1e-12,
      TRUE ~ observed == expected
    )
  )

if (!all(metadata_validation$matches)) {
  print(metadata_validation, n = Inf, width = Inf)
  stop("Metadata validation failed.")
}

readr::write_csv(
  metadata_validation,
  file.path(tables_dir, "metadata_validation_sex_reframe.csv")
)

metadata_with_cluster <- metadata |>
  dplyr::select(
    participant_id,
    age,
    sex,
    pubertal_development,
    endocrine_status,
    uctd_status,
    periodontal_status,
    baseline_cluster
  )

readr::write_csv(
  metadata_with_cluster,
  file.path(tables_dir, "pediatric_metadata_with_baseline_cluster.csv")
)

cohort_metadata_summary <- dplyr::bind_rows(
  tibble::tibble(
    variable = "age",
    category = "overall",
    n = nrow(metadata),
    mean = mean(metadata$age),
    sd = stats::sd(metadata$age),
    median = stats::median(metadata$age),
    min = min(metadata$age),
    max = max(metadata$age)
  ),
  categorical_overall_summary(metadata, "sex"),
  categorical_overall_summary(metadata, "pubertal_development"),
  categorical_overall_summary(metadata, "endocrine_status"),
  categorical_overall_summary(metadata, "uctd_status"),
  categorical_overall_summary(metadata, "periodontal_status")
)

readr::write_csv(
  cohort_metadata_summary,
  file.path(tables_dir, "cohort_metadata_summary.csv")
)


# ------------------------------------------------------------
# 4. Sex association with the fixed baseline partition
# ------------------------------------------------------------

sex_cluster_tab <- table(metadata$sex, metadata$baseline_cluster)
sex_cluster_fisher <- stats::fisher.test(sex_cluster_tab)
sex_cluster_v <- cramers_v(metadata$sex, metadata$baseline_cluster)
sex_cluster_v_ci <- bootstrap_cramers_v(
  x = metadata$sex,
  y = metadata$baseline_cluster,
  B = N_BOOT_EFFECT,
  seed = SEED + 1L
)

sex_cluster_counts <- as.data.frame(sex_cluster_tab) |>
  tibble::as_tibble() |>
  dplyr::rename(
    sex = Var1,
    baseline_cluster = Var2,
    n = Freq
  ) |>
  dplyr::mutate(
    sex = as.character(sex),
    baseline_cluster = as.integer(as.character(baseline_cluster))
  )

sex_cluster_association <- tibble::tibble(
  n = sum(sex_cluster_tab),
  fisher_p_value = sex_cluster_fisher$p.value,
  cramers_v = sex_cluster_v,
  cramers_v_ci_lower = sex_cluster_v_ci[["lower"]],
  cramers_v_ci_upper = sex_cluster_v_ci[["upper"]],
  valid_bootstrap = as.integer(sex_cluster_v_ci[["valid_bootstrap"]])
)

readr::write_csv(
  sex_cluster_counts,
  file.path(tables_dir, "sex_by_baseline_cluster_counts.csv")
)

readr::write_csv(
  sex_cluster_association,
  file.path(tables_dir, "sex_baseline_cluster_association.csv")
)


# ------------------------------------------------------------
# 5. Sex associations with PC1-PC3, including BH correction
# ------------------------------------------------------------

scores_df <- tibble::tibble(
  participant_id = sample_ids,
  PC1 = scores_baseline[, 1],
  PC2 = scores_baseline[, 2],
  PC3 = scores_baseline[, 3]
) |>
  dplyr::left_join(
    metadata |>
      dplyr::select(participant_id, sex, baseline_cluster),
    by = "participant_id"
  )

pc_names <- c("PC1", "PC2", "PC3")

pc_sex_results <- lapply(
  seq_along(pc_names),
  function(i) {
    pc_name <- pc_names[[i]]
    x <- scores_df[[pc_name]]
    sex <- scores_df$sex

    wilcox_res <- stats::wilcox.test(
      x ~ sex,
      conf.int = TRUE
    )

    x_f <- x[sex == "F"]
    x_m <- x[sex == "M"]

    rank_biserial <- rank_biserial_two_group(x_f, x_m)
    rank_biserial_ci <- bootstrap_two_group_effect(
      x = x,
      sex = sex,
      effect_fun = rank_biserial_two_group,
      B = N_BOOT_EFFECT,
      seed = SEED + 100L + i
    )

    hedges_g <- hedges_g_two_group(x_f, x_m)
    hedges_g_ci <- bootstrap_two_group_effect(
      x = x,
      sex = sex,
      effect_fun = hedges_g_two_group,
      B = N_BOOT_EFFECT,
      seed = SEED + 200L + i
    )

    tibble::tibble(
      PC = pc_name,
      n_F = length(x_f),
      n_M = length(x_m),
      mean_F = mean(x_f),
      mean_M = mean(x_m),
      median_F = stats::median(x_f),
      median_M = stats::median(x_m),
      wilcoxon_W = as.numeric(wilcox_res$statistic),
      p_value_raw = wilcox_res$p.value,
      hodges_lehmann_F_minus_M = as.numeric(wilcox_res$estimate),
      hodges_lehmann_ci_lower = as.numeric(wilcox_res$conf.int[[1]]),
      hodges_lehmann_ci_upper = as.numeric(wilcox_res$conf.int[[2]]),
      rank_biserial_F_vs_M = rank_biserial,
      rank_biserial_ci_lower = rank_biserial_ci[["lower"]],
      rank_biserial_ci_upper = rank_biserial_ci[["upper"]],
      hedges_g_F_minus_M = hedges_g,
      hedges_g_ci_lower = hedges_g_ci[["lower"]],
      hedges_g_ci_upper = hedges_g_ci[["upper"]]
    )
  }
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    p_value_BH = stats::p.adjust(p_value_raw, method = "BH"),
    significant_BH_0_05 = p_value_BH < 0.05
  )

readr::write_csv(
  pc_sex_results,
  file.path(tables_dir, "PC1_PC3_sex_associations_BH.csv")
)

p_pc1_sex <- ggplot2::ggplot(
  scores_df,
  ggplot2::aes(x = sex, y = PC1)
) +
  ggplot2::geom_boxplot(width = 0.5, outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.08, height = 0, size = 2) +
  ggplot2::labs(x = "Sex", y = "PC1 score") +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "diagnostic_PC1_by_sex.png"),
  p_pc1_sex,
  width = 6,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 6. Feature-wise sex residualization
# ------------------------------------------------------------

# The matrix formulation below is algebraically equivalent to
# fitting 793 separate models with the same design matrix.
design_sex <- stats::model.matrix(~ sex, data = metadata)
coef_sex <- qr.solve(design_sex, X_snv)
X_sex_residualized <- X_snv - design_sex %*% coef_sex

rownames(X_sex_residualized) <- sample_ids
colnames(X_sex_residualized) <- colnames(X_snv)

max_abs_residual_orthogonality <- max(
  abs(crossprod(design_sex, X_sex_residualized))
)

pca_sex_residualized <- FactoMineR::PCA(
  X_sex_residualized,
  scale.unit = FALSE,
  ncp = N_PC_BASELINE,
  graph = FALSE
)

scores_sex_residualized <- as.matrix(
  pca_sex_residualized$ind$coord[
    ,
    seq_len(N_PC_BASELINE),
    drop = FALSE
  ]
)
rownames(scores_sex_residualized) <- sample_ids

cluster_sex_residualized_raw <- ward_clusters(
  scores_sex_residualized,
  k = K_BASELINE
)
names(cluster_sex_residualized_raw) <- sample_ids

cluster_sex_residualized <- align_cluster_labels(
  reference = cluster_baseline,
  candidate = cluster_sex_residualized_raw
)

sex_residualized_ARI <- mclust::adjustedRandIndex(
  cluster_baseline,
  cluster_sex_residualized_raw
)

sex_residualized_silhouette <- cluster::silhouette(
  cluster_sex_residualized_raw,
  stats::dist(scores_sex_residualized, method = "euclidean")
)

sex_residualized_average_silhouette <- mean(
  sex_residualized_silhouette[, "sil_width"]
)

sex_residualized_negative_silhouette_n <- sum(
  sex_residualized_silhouette[, "sil_width"] < 0
)

resid_eig <- as.data.frame(pca_sex_residualized$eig)
resid_cumulative_variance_3PC <- resid_eig[N_PC_BASELINE, 3]


# ------------------------------------------------------------
# 7. Diagnostic sex association after residualization
# ------------------------------------------------------------

sex_residualized_factor <- factor(
  cluster_sex_residualized[metadata$participant_id],
  levels = 1:K_BASELINE
)

sex_residualized_tab <- table(metadata$sex, sex_residualized_factor)
sex_residualized_fisher <- stats::fisher.test(sex_residualized_tab)
sex_residualized_v <- cramers_v(metadata$sex, sex_residualized_factor)
sex_residualized_v_ci <- bootstrap_cramers_v(
  x = metadata$sex,
  y = sex_residualized_factor,
  B = N_BOOT_EFFECT,
  seed = SEED + 300L
)

sex_residualized_assignments <- tibble::tibble(
  participant_id = sample_ids,
  sex = as.character(metadata$sex),
  baseline_cluster = as.integer(cluster_baseline[sample_ids]),
  sex_residualized_cluster = as.integer(
    cluster_sex_residualized[sample_ids]
  ),
  assignment_changed = baseline_cluster != sex_residualized_cluster,
  residualized_PC1 = scores_sex_residualized[, 1],
  residualized_PC2 = scores_sex_residualized[, 2],
  residualized_PC3 = scores_sex_residualized[, 3]
)

sex_residualized_cluster_sizes <- as.integer(
  table(
    factor(
      cluster_sex_residualized,
      levels = 1:K_BASELINE
    )
  )
)

sex_residualized_summary <- tibble::tibble(
  n_participants = nrow(X_snv),
  baseline_cluster_sizes = paste(
    as.integer(table(cluster_baseline)),
    collapse = "/"
  ),
  residualized_cluster_sizes = paste(
    sex_residualized_cluster_sizes,
    collapse = "/"
  ),
  ARI_vs_baseline = sex_residualized_ARI,
  participants_changing_aligned_assignment = sum(
    sex_residualized_assignments$assignment_changed
  ),
  residualized_average_silhouette =
    sex_residualized_average_silhouette,
  residualized_negative_silhouette_n =
    sex_residualized_negative_silhouette_n,
  residualized_3PC_cumulative_variance_percent =
    resid_cumulative_variance_3PC,
  sex_cluster_fisher_p_after_residualization =
    sex_residualized_fisher$p.value,
  sex_cluster_cramers_v_after_residualization =
    sex_residualized_v,
  sex_cluster_cramers_v_ci_lower =
    sex_residualized_v_ci[["lower"]],
  sex_cluster_cramers_v_ci_upper =
    sex_residualized_v_ci[["upper"]],
  max_abs_design_residual_crossproduct =
    max_abs_residual_orthogonality
)

baseline_vs_sex_residualized <- as.data.frame.matrix(
  table(
    Baseline = sex_residualized_assignments$baseline_cluster,
    SexResidualized =
      sex_residualized_assignments$sex_residualized_cluster
  )
) |>
  tibble::rownames_to_column(var = "baseline_cluster")

readr::write_csv(
  sex_residualized_assignments,
  file.path(tables_dir, "sex_residualized_cluster_assignments.csv")
)

readr::write_csv(
  sex_residualized_summary,
  file.path(tables_dir, "sex_residualized_clustering_summary.csv")
)

readr::write_csv(
  baseline_vs_sex_residualized,
  file.path(tables_dir, "baseline_vs_sex_residualized_clusters.csv")
)

saveRDS(
  X_sex_residualized,
  file.path(processed_dir, "pediatric_X_snv_sex_residualized.rds")
)

saveRDS(
  pca_sex_residualized,
  file.path(processed_dir, "pediatric_PCA_sex_residualized.rds")
)

saveRDS(
  cluster_sex_residualized,
  file.path(processed_dir, "pediatric_clusters_sex_residualized.rds")
)


# ------------------------------------------------------------
# 8. Compact reporting table and result bundle
# ------------------------------------------------------------

sex_sensitivity_summary <- dplyr::bind_rows(
  tibble::tibble(
    section = "Baseline sex association",
    metric = c(
      "Sex x baseline cluster Fisher p",
      "Sex x baseline cluster Cramer's V",
      "Sex x baseline cluster Cramer's V CI lower",
      "Sex x baseline cluster Cramer's V CI upper"
    ),
    value_numeric = c(
      sex_cluster_association$fisher_p_value,
      sex_cluster_association$cramers_v,
      sex_cluster_association$cramers_v_ci_lower,
      sex_cluster_association$cramers_v_ci_upper
    )
  ),
  pc_sex_results |>
    dplyr::transmute(
      section = paste0(PC, " sex association"),
      metric = "Wilcoxon raw p",
      value_numeric = p_value_raw
    ),
  pc_sex_results |>
    dplyr::transmute(
      section = paste0(PC, " sex association"),
      metric = "Wilcoxon BH-adjusted p",
      value_numeric = p_value_BH
    ),
  pc_sex_results |>
    dplyr::transmute(
      section = paste0(PC, " sex association"),
      metric = "Rank-biserial correlation F vs M",
      value_numeric = rank_biserial_F_vs_M
    ),
  pc_sex_results |>
    dplyr::transmute(
      section = paste0(PC, " sex association"),
      metric = "Hedges g F minus M",
      value_numeric = hedges_g_F_minus_M
    ),
  tibble::tibble(
    section = "Sex-residualized sensitivity",
    metric = c(
      "ARI vs baseline",
      "Participants changing aligned assignment",
      "Residualized average silhouette",
      "Residualized 3-PC cumulative variance (%)",
      "Sex x residualized cluster Fisher p",
      "Sex x residualized cluster Cramer's V",
      "Sex x residualized cluster Cramer's V CI lower",
      "Sex x residualized cluster Cramer's V CI upper"
    ),
    value_numeric = c(
      sex_residualized_summary$ARI_vs_baseline,
      sex_residualized_summary$participants_changing_aligned_assignment,
      sex_residualized_summary$residualized_average_silhouette,
      sex_residualized_summary$residualized_3PC_cumulative_variance_percent,
      sex_residualized_summary$sex_cluster_fisher_p_after_residualization,
      sex_residualized_summary$sex_cluster_cramers_v_after_residualization,
      sex_residualized_summary$sex_cluster_cramers_v_ci_lower,
      sex_residualized_summary$sex_cluster_cramers_v_ci_upper
    )
  )
) |>
  dplyr::mutate(value_text = format(value_numeric, digits = 10))

readr::write_csv(
  sex_sensitivity_summary,
  file.path(tables_dir, "sex_sensitivity_summary.csv")
)

sex_sensitivity_objects <- list(
  settings = list(
    seed = SEED,
    baseline_n_pc = N_PC_BASELINE,
    baseline_k = K_BASELINE,
    n_boot_effect = N_BOOT_EFFECT,
    pc_multiplicity_adjustment = "Benjamini-Hochberg across PC1-PC3"
  ),
  metadata = metadata_with_cluster,
  metadata_validation = metadata_validation,
  cohort_metadata_summary = cohort_metadata_summary,
  sex_cluster_counts = sex_cluster_counts,
  sex_cluster_association = sex_cluster_association,
  pc_sex_results = pc_sex_results,
  sex_residualized_assignments = sex_residualized_assignments,
  sex_residualized_summary = sex_residualized_summary,
  baseline_vs_sex_residualized = baseline_vs_sex_residualized,
  sex_sensitivity_summary = sex_sensitivity_summary
)

saveRDS(
  sex_sensitivity_objects,
  file.path(processed_dir, "sex_sensitivity_objects.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(tables_dir, "sessionInfo_sex_sensitivity.txt")
)


# ------------------------------------------------------------
# 9. Console report
# ------------------------------------------------------------

pc1_row <- pc_sex_results |>
  dplyr::filter(PC == "PC1")

cat("\n========================================\n")
cat("SEX-SENSITIVITY ANALYSIS COMPLETE\n")
cat("========================================\n")
cat(
  "Sex x baseline cluster Fisher p: ",
  signif(sex_cluster_association$fisher_p_value, 4),
  "\n",
  sep = ""
)
cat(
  "PC1 sex Wilcoxon raw p / BH p: ",
  signif(pc1_row$p_value_raw, 4),
  " / ",
  signif(pc1_row$p_value_BH, 4),
  "\n",
  sep = ""
)
cat(
  "Sex-residualized ARI vs baseline: ",
  round(sex_residualized_ARI, 3),
  "\n",
  sep = ""
)
cat(
  "Sex-residualized aligned cluster sizes: ",
  sex_residualized_summary$residualized_cluster_sizes,
  "\n",
  sep = ""
)
cat(
  "Participants changing aligned assignment: ",
  sex_residualized_summary$participants_changing_aligned_assignment,
  "\n",
  sep = ""
)
cat("========================================\n")
