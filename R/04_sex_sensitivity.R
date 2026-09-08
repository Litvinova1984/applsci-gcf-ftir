# ============================================================
# 04_sex_sensitivity.R
#
# Pediatric GCF FTIR
# Clinical-metadata associations and sex-sensitivity analysis
#
# PURPOSE
# -------
# This script evaluates how the baseline three-cluster solution
# relates to the available participant-level clinical metadata,
# with particular attention to sex because sex is associated
# with the dominant spectral axis (PC1) in this cohort.
#
# The script addresses three questions:
#
#   1. Are the baseline spectral clusters associated with the
#      available clinical metadata?
#        - sex, endocrine status, UCTD status, and periodontal
#          status: Fisher's exact tests
#        - age: Kruskal-Wallis test
#
#   2. Do the first three PCA scores differ by sex?
#        - Wilcoxon rank-sum tests for PC1-PC3
#        - effect sizes and bootstrap 95% confidence intervals
#
#   3. How dependent is the baseline cluster partition on
#      sex-associated spectral variation?
#        - remove the linear sex-associated mean component from
#          each SNV-normalized spectral variable
#        - recompute PCA using the same three-PC specification
#        - recompute Ward.D2 clustering at k = 3
#        - compare the residualized and baseline partitions
#          using the adjusted Rand index (ARI)
#
# IMPORTANT INTERPRETIVE POINT
# ----------------------------
# The residualization analysis is a sensitivity analysis, not
# a causal model. It asks how much the observed cluster
# structure changes after removing spectral variation that is
# linearly associated with sex in this sample. A change in the
# partition indicates that sex-associated variation contributes
# to the observed spectral organization; it does not establish
# a biological mechanism or imply that the clusters are simply
# equivalent to sex categories.
#
# Pubertal development is retained in the public metadata table
# and summarized descriptively. It is not tested inferentially
# because 13 of 25 participants have unknown pubertal status,
# leaving sparse observed categories in this small cohort.
#
# All random resampling uses the project-wide SEED defined in
# R/00_utils.R. Confidence intervals are percentile bootstrap
# intervals and should be interpreted cautiously given n = 25.
# ============================================================


# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}

source(
  file.path(
    "R",
    "00_utils.R"
  )
)


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
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 1. Paths and fixed analytical settings
# ------------------------------------------------------------

processed_dir <- file.path(
  "data",
  "pediatric",
  "processed"
)

metadata_file <- file.path(
  "data",
  "pediatric",
  "pediatric_metadata.csv"
)

x_file <- file.path(
  processed_dir,
  "pediatric_X_snv.rds"
)

baseline_cluster_file <- file.path(
  processed_dir,
  "pediatric_clusters_baseline.rds"
)

tables_dir <- "tables"
figures_dir <- "figures"

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


if (!file.exists(metadata_file)) {
  stop(
    "pediatric_metadata.csv not found at ",
    metadata_file
  )
}

if (!file.exists(x_file)) {
  stop(
    "pediatric_X_snv.rds not found. ",
    "Run R/01_validate_public_inputs.R and verify the public processed inputs."
  )
}

if (!file.exists(baseline_cluster_file)) {
  stop(
    "pediatric_clusters_baseline.rds not found. ",
    "Run R/02_pca_clustering.R first."
  )
}


N_PC_BASELINE <- 3L
K_BASELINE <- 3L

# Number of bootstrap replicates used for effect-size confidence
# intervals. This affects the precision of the percentile CIs, not
# the point estimates or hypothesis-test p-values.
N_BOOT_EFFECT <- 5000L


# ------------------------------------------------------------
# 2. Load and validate data
# ------------------------------------------------------------

metadata <- readr::read_csv(
  metadata_file,
  show_col_types = FALSE
) |>
  dplyr::mutate(
    participant_id =
      as.character(participant_id),
    age =
      as.numeric(age),
    sex =
      factor(
        sex,
        levels = c("F", "M")
      ),
    pubertal_development =
      factor(
        pubertal_development,
        levels = c(
          "Normal",
          "Early",
          "Delayed",
          "Unknown"
        )
      ),
    endocrine_status =
      factor(
        endocrine_status,
        levels = c(
          "No documented pathology",
          "Documented pathology",
          "Unknown"
        )
      ),
    uctd_status =
      factor(
        uctd_status,
        levels = c(
          "No documented pathology",
          "Documented pathology"
        )
      ),
    periodontal_status =
      factor(
        periodontal_status,
        levels = c(
          "No documented pathology",
          "Documented pathology"
        )
      )
  )


X_snv <- as.matrix(
  readRDS(x_file)
)

sample_ids <- rownames(X_snv)

if (is.null(sample_ids)) {
  stop("X_snv must have participant IDs as row names.")
}

if (!identical(
  dim(X_snv),
  c(25L, 793L)
)) {
  stop(
    "Expected X_snv dimensions 25 x 793; obtained ",
    paste(dim(X_snv), collapse = " x "),
    "."
  )
}

if (any(!is.finite(X_snv))) {
  stop("Non-finite values detected in X_snv.")
}

if (nrow(metadata) != 25L) {
  stop(
    "Expected 25 metadata rows; obtained ",
    nrow(metadata),
    "."
  )
}

if (anyDuplicated(metadata$participant_id) > 0) {
  stop("Duplicate participant IDs found in metadata.")
}

if (!setequal(
  metadata$participant_id,
  sample_ids
)) {
  stop(
    "Participant IDs in pediatric_metadata.csv ",
    "do not match X_snv row names."
  )
}


# Reorder metadata to exactly match the spectral matrix.

metadata <- metadata[
  match(
    sample_ids,
    metadata$participant_id
  ),
  ,
  drop = FALSE
]


cluster_baseline_raw <- readRDS(
  baseline_cluster_file
)

if (length(cluster_baseline_raw) != nrow(X_snv)) {
  stop(
    "Baseline cluster vector length does not match X_snv."
  )
}

if (is.null(names(cluster_baseline_raw))) {

  cluster_baseline <- as.integer(
    cluster_baseline_raw
  )

  names(cluster_baseline) <- sample_ids

} else {

  if (!setequal(
    names(cluster_baseline_raw),
    sample_ids
  )) {
    stop(
      "Participant IDs in the baseline cluster vector ",
      "do not match X_snv."
    )
  }

  cluster_baseline_raw <-
    cluster_baseline_raw[
      sample_ids
    ]

  if (is.factor(cluster_baseline_raw)) {
    cluster_baseline <- as.integer(
      as.character(
        cluster_baseline_raw
      )
    )
  } else {
    cluster_baseline <- as.integer(
      cluster_baseline_raw
    )
  }

  names(cluster_baseline) <- sample_ids
}


if (length(unique(cluster_baseline)) != K_BASELINE) {
  stop(
    "Expected a three-cluster baseline solution."
  )
}


metadata$cluster <- factor(
  cluster_baseline[
    metadata$participant_id
  ],
  levels = 1:K_BASELINE
)


# ------------------------------------------------------------
# 3. Analysis-cohort integrity checks
# ------------------------------------------------------------

# These fixed values define the frozen 25-participant analysis
# cohort distributed with the reproducibility repository. They act
# as integrity checks against accidental row loss, recoding, or use
# of a different metadata file.
metadata_validation <- tibble::tibble(
  check = c(
    "n_participants",
    "age_min",
    "age_max",
    "age_mean",
    "age_sd",
    "sex_F",
    "sex_M",
    "pubertal_Normal",
    "pubertal_Early",
    "pubertal_Delayed",
    "pubertal_Unknown",
    "endocrine_documented",
    "endocrine_no_documented",
    "endocrine_unknown",
    "uctd_documented",
    "uctd_no_documented",
    "periodontal_documented",
    "periodontal_no_documented"
  ),
  observed = c(
    nrow(metadata),
    min(metadata$age),
    max(metadata$age),
    mean(metadata$age),
    stats::sd(metadata$age),
    sum(metadata$sex == "F"),
    sum(metadata$sex == "M"),
    sum(metadata$pubertal_development == "Normal"),
    sum(metadata$pubertal_development == "Early"),
    sum(metadata$pubertal_development == "Delayed"),
    sum(metadata$pubertal_development == "Unknown"),
    sum(metadata$endocrine_status == "Documented pathology"),
    sum(metadata$endocrine_status == "No documented pathology"),
    sum(metadata$endocrine_status == "Unknown"),
    sum(metadata$uctd_status == "Documented pathology"),
    sum(metadata$uctd_status == "No documented pathology"),
    sum(metadata$periodontal_status == "Documented pathology"),
    sum(metadata$periodontal_status == "No documented pathology")
  ),
  expected = c(
    25,
    6,
    17,
    11.92,
    3.37,
    12,
    13,
    8,
    3,
    1,
    13,
    19,
    5,
    1,
    22,
    3,
    3,
    22
  )
) |>
  dplyr::mutate(
    match = dplyr::case_when(
      check %in% c(
        "age_mean",
        "age_sd"
      ) ~
        abs(
          round(observed, 2) -
            expected
        ) < 1e-12,
      TRUE ~
        observed == expected
    )
  )


if (!all(metadata_validation$match)) {
  print(
    metadata_validation,
    n = Inf,
    width = Inf
  )
  stop(
    "Metadata validation failed."
  )
}


cat("\n========================================\n")
cat("METADATA VALIDATION PASSED\n")
cat("========================================\n")

print(
  metadata_validation,
  n = Inf,
  width = Inf
)


readr::write_csv(
  metadata_validation,
  file.path(
    tables_dir,
    "metadata_validation.csv"
  )
)


# Save a public analysis-ready table including only
# anonymized analytical variables and the baseline cluster.

metadata_with_cluster <- metadata |>
  dplyr::select(
    participant_id,
    age,
    sex,
    pubertal_development,
    endocrine_status,
    uctd_status,
    periodontal_status,
    cluster
  )

readr::write_csv(
  metadata_with_cluster,
  file.path(
    tables_dir,
    "pediatric_metadata_with_cluster.csv"
  )
)


# ------------------------------------------------------------
# 4. Helper functions for effect sizes and confidence intervals
# ------------------------------------------------------------

cramers_v <- function(x, y) {

  tab <- table(
    x,
    y
  )

  tab <- tab[
    rowSums(tab) > 0,
    colSums(tab) > 0,
    drop = FALSE
  ]

  if (
    nrow(tab) < 2 ||
      ncol(tab) < 2
  ) {
    return(NA_real_)
  }

  n <- sum(tab)

  expected <- outer(
    rowSums(tab),
    colSums(tab)
  ) / n

  valid <- expected > 0

  chi_sq <- sum(
    (
      tab[valid] -
        expected[valid]
    )^2 /
      expected[valid]
  )

  denom <- min(
    nrow(tab) - 1,
    ncol(tab) - 1
  )

  sqrt(
    chi_sq /
      (
        n *
          denom
      )
  )
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

      idx <- sample.int(
        n,
        size = n,
        replace = TRUE
      )

      cramers_v(
        x[idx],
        y[idx]
      )
    }
  )

  vals <- vals[
    is.finite(vals)
  ]

  c(
    lower =
      as.numeric(
        stats::quantile(
          vals,
          0.025,
          na.rm = TRUE
        )
      ),
    upper =
      as.numeric(
        stats::quantile(
          vals,
          0.975,
          na.rm = TRUE
        )
      ),
    valid_bootstrap =
      length(vals)
  )
}


epsilon_squared_kw <- function(
    x,
    group
) {

  keep <- is.finite(x) &
    !is.na(group)

  x <- x[keep]
  group <- droplevels(
    factor(group[keep])
  )

  if (nlevels(group) < 2) {
    return(NA_real_)
  }

  kw <- stats::kruskal.test(
    x ~ group
  )

  H <- as.numeric(
    kw$statistic
  )

  n <- length(x)
  k <- nlevels(group)

  eps2 <- (
    H - k + 1
  ) /
    (
      n - k
    )

  max(
    0,
    eps2
  )
}


bootstrap_epsilon_squared <- function(
    x,
    group,
    B = N_BOOT_EFFECT,
    seed = SEED
) {

  n <- length(x)

  set.seed(seed)

  vals <- replicate(
    B,
    {

      idx <- sample.int(
        n,
        size = n,
        replace = TRUE
      )

      epsilon_squared_kw(
        x[idx],
        group[idx]
      )
    }
  )

  vals <- vals[
    is.finite(vals)
  ]

  c(
    lower =
      as.numeric(
        stats::quantile(
          vals,
          0.025,
          na.rm = TRUE
        )
      ),
    upper =
      as.numeric(
        stats::quantile(
          vals,
          0.975,
          na.rm = TRUE
        )
      ),
    valid_bootstrap =
      length(vals)
  )
}


rank_biserial_two_group <- function(
    x_f,
    x_m
) {

  x_f <- x_f[
    is.finite(x_f)
  ]

  x_m <- x_m[
    is.finite(x_m)
  ]

  n_f <- length(x_f)
  n_m <- length(x_m)

  pooled <- c(
    x_f,
    x_m
  )

  ranks <- rank(
    pooled,
    ties.method = "average"
  )

  R_f <- sum(
    ranks[
      seq_len(n_f)
    ]
  )

  U_f <- R_f -
    n_f *
    (
      n_f + 1
    ) / 2

  2 *
    U_f /
    (
      n_f *
        n_m
    ) -
    1
}


hedges_g_two_group <- function(
    x_f,
    x_m
) {

  x_f <- x_f[
    is.finite(x_f)
  ]

  x_m <- x_m[
    is.finite(x_m)
  ]

  n_f <- length(x_f)
  n_m <- length(x_m)

  df <- n_f + n_m - 2

  pooled_sd <- sqrt(
    (
      (n_f - 1) *
        stats::var(x_f) +
        (n_m - 1) *
        stats::var(x_m)
    ) /
      df
  )

  if (
    !is.finite(pooled_sd) ||
      pooled_sd == 0
  ) {
    return(NA_real_)
  }

  d <- (
    mean(x_f) -
      mean(x_m)
  ) /
    pooled_sd

  J <- 1 -
    3 /
    (
      4 *
        (
          n_f + n_m
        ) -
        9
    )

  J * d
}


bootstrap_two_group_effect <- function(
    x,
    sex,
    effect_fun,
    B = N_BOOT_EFFECT,
    seed = SEED
) {

  x_f <- x[
    sex == "F"
  ]

  x_m <- x[
    sex == "M"
  ]

  set.seed(seed)

  vals <- replicate(
    B,
    {

      boot_f <- sample(
        x_f,
        size = length(x_f),
        replace = TRUE
      )

      boot_m <- sample(
        x_m,
        size = length(x_m),
        replace = TRUE
      )

      effect_fun(
        boot_f,
        boot_m
      )
    }
  )

  vals <- vals[
    is.finite(vals)
  ]

  c(
    lower =
      as.numeric(
        stats::quantile(
          vals,
          0.025,
          na.rm = TRUE
        )
      ),
    upper =
      as.numeric(
        stats::quantile(
          vals,
          0.975,
          na.rm = TRUE
        )
      ),
    valid_bootstrap =
      length(vals)
  )
}


ward_clusters <- function(
    scores,
    k
) {

  hc <- stats::hclust(
    stats::dist(
      scores,
      method = "euclidean"
    ),
    method = "ward.D2"
  )

  stats::cutree(
    hc,
    k = k
  )
}


all_permutations <- function(x) {

  if (length(x) == 1L) {
    return(
      list(x)
    )
  }

  out <- list()

  counter <- 1L

  for (i in seq_along(x)) {

    sub_perm <- all_permutations(
      x[-i]
    )

    for (p in sub_perm) {

      out[[counter]] <- c(
        x[i],
        p
      )

      counter <- counter + 1L
    }
  }

  out
}


align_cluster_labels <- function(
    reference,
    candidate
) {

  ref_levels <- sort(
    unique(reference)
  )

  cand_levels <- sort(
    unique(candidate)
  )

  if (
    length(ref_levels) !=
      length(cand_levels)
  ) {
    stop(
      "Reference and candidate partitions ",
      "have different numbers of clusters."
    )
  }

  permutations <- all_permutations(
    ref_levels
  )

  agreements <- vapply(
    permutations,
    function(p) {

      map <- stats::setNames(
        p,
        cand_levels
      )

      aligned <- unname(
        map[
          as.character(
            candidate
          )
        ]
      )

      sum(
        aligned == reference
      )
    },
    numeric(1)
  )

  best <- which.max(
    agreements
  )

  map <- stats::setNames(
    permutations[[best]],
    cand_levels
  )

  aligned <- as.integer(
    unname(
      map[
        as.character(
          candidate
        )
      ]
    )
  )

  names(aligned) <- names(candidate)

  aligned
}


# ------------------------------------------------------------
# Main outputs created by this script
# ------------------------------------------------------------
#
# tables/
#   clinical_age_by_cluster.csv
#   clinical_categorical_counts_by_cluster.csv
#   clinical_categorical_tests_effect_sizes.csv
#   clinical_age_test_effect_size.csv
#   clinical_PC1_PC3_sex_effect_sizes.csv
#   sex_residualized_cluster_assignments.csv
#   sex_residualized_clustering_summary.csv
#   baseline_vs_sex_residualized_clusters.csv
#   sex_sensitivity_summary.csv
#   sessionInfo_sex_sensitivity.txt
#
# data/pediatric/processed/
#   pediatric_X_snv_sex_residualized.rds
#   pediatric_PCA_sex_residualized.rds
#   pediatric_clusters_sex_residualized.rds
#   sex_sensitivity_objects.rds
#
# figures/
#   clinical_PC1_by_sex.png
#
# ------------------------------------------------------------


# ============================================================
# PART A. CLINICAL-METADATA ASSOCIATIONS WITH BASELINE CLUSTERS
# ============================================================


# ------------------------------------------------------------
# 5. Descriptive metadata by baseline cluster
# ------------------------------------------------------------

age_by_cluster <- metadata |>
  dplyr::group_by(
    cluster
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_age =
      mean(age),
    sd_age =
      stats::sd(age),
    median_age =
      stats::median(age),
    q25_age =
      as.numeric(
        stats::quantile(
          age,
          0.25
        )
      ),
    q75_age =
      as.numeric(
        stats::quantile(
          age,
          0.75
        )
      ),
    min_age =
      min(age),
    max_age =
      max(age),
    .groups = "drop"
  )


categorical_cluster_counts <- dplyr::bind_rows(
  metadata |>
    dplyr::count(
      cluster,
      sex,
      name = "n"
    ) |>
    dplyr::transmute(
      variable = "sex",
      cluster,
      category =
        as.character(sex),
      n
    ),
  metadata |>
    dplyr::count(
      cluster,
      pubertal_development,
      name = "n"
    ) |>
    dplyr::transmute(
      variable =
        "pubertal_development",
      cluster,
      category =
        as.character(
          pubertal_development
        ),
      n
    ),
  metadata |>
    dplyr::count(
      cluster,
      endocrine_status,
      name = "n"
    ) |>
    dplyr::transmute(
      variable =
        "endocrine_status",
      cluster,
      category =
        as.character(
          endocrine_status
        ),
      n
    ),
  metadata |>
    dplyr::count(
      cluster,
      uctd_status,
      name = "n"
    ) |>
    dplyr::transmute(
      variable =
        "uctd_status",
      cluster,
      category =
        as.character(
          uctd_status
        ),
      n
    ),
  metadata |>
    dplyr::count(
      cluster,
      periodontal_status,
      name = "n"
    ) |>
    dplyr::transmute(
      variable =
        "periodontal_status",
      cluster,
      category =
        as.character(
          periodontal_status
        ),
      n
    )
)


readr::write_csv(
  age_by_cluster,
  file.path(
    tables_dir,
    "clinical_age_by_cluster.csv"
  )
)

readr::write_csv(
  categorical_cluster_counts,
  file.path(
    tables_dir,
    "clinical_categorical_counts_by_cluster.csv"
  )
)


# ------------------------------------------------------------
# 6. Fisher exact tests for categorical metadata
#
# Fisher's exact test is used because several contingency-table
# cells are small. For endocrine status, the single "Unknown"
# participant is retained as an explicit analytical category rather
# than being silently excluded.
#
# Cramer's V is reported as an effect-size measure:
#   0 indicates no association; larger values indicate a stronger
#   association. Bootstrap CIs quantify uncertainty in the effect
#   size estimate in this small cohort.
# ------------------------------------------------------------

categorical_vars <- c(
  "sex",
  "endocrine_status",
  "uctd_status",
  "periodontal_status"
)

categorical_tests <- lapply(
  seq_along(categorical_vars),
  function(i) {

    var_name <-
      categorical_vars[i]

    x <- metadata[[var_name]]

    tab <- table(
      x,
      metadata$cluster
    )

    fisher_res <- stats::fisher.test(
      tab
    )

    v <- cramers_v(
      x,
      metadata$cluster
    )

    v_ci <- bootstrap_cramers_v(
      x = x,
      y = metadata$cluster,
      B = N_BOOT_EFFECT,
      seed = SEED + i
    )

    tibble::tibble(
      variable = var_name,
      test = "Fisher exact",
      n = sum(tab),
      table_dimensions =
        paste(
          dim(tab),
          collapse = "x"
        ),
      p_value =
        fisher_res$p.value,
      cramers_v = v,
      cramers_v_ci_lower =
        v_ci[
          "lower"
        ],
      cramers_v_ci_upper =
        v_ci[
          "upper"
        ],
      valid_bootstrap =
        as.integer(
          v_ci[
            "valid_bootstrap"
          ]
        )
    )
  }
)

categorical_tests <- dplyr::bind_rows(
  categorical_tests
)


# ------------------------------------------------------------
# 7. Kruskal-Wallis test for age
#
# Age is treated as a quantitative variable. The Kruskal-Wallis
# test compares age distributions across the three clusters without
# assuming normality. Epsilon-squared is reported as a nonparametric
# effect-size estimate; small negative finite-sample estimates are
# truncated to zero.
# ------------------------------------------------------------

age_kw <- stats::kruskal.test(
  age ~ cluster,
  data = metadata
)

age_eps2 <- epsilon_squared_kw(
  metadata$age,
  metadata$cluster
)

age_eps2_ci <- bootstrap_epsilon_squared(
  x = metadata$age,
  group = metadata$cluster,
  B = N_BOOT_EFFECT,
  seed = SEED + 100L
)


age_test <- tibble::tibble(
  variable = "age",
  test = "Kruskal-Wallis",
  n = nrow(metadata),
  statistic =
    as.numeric(
      age_kw$statistic
    ),
  df =
    as.numeric(
      age_kw$parameter
    ),
  p_value =
    age_kw$p.value,
  epsilon_squared =
    age_eps2,
  epsilon_squared_ci_lower =
    age_eps2_ci[
      "lower"
    ],
  epsilon_squared_ci_upper =
    age_eps2_ci[
      "upper"
    ],
  valid_bootstrap =
    as.integer(
      age_eps2_ci[
        "valid_bootstrap"
      ]
    )
)


cat("\n========================================\n")
cat("CLINICAL METADATA ASSOCIATIONS\n")
cat("========================================\n")

print(
  categorical_tests,
  n = Inf,
  width = Inf
)

print(
  age_test,
  n = Inf,
  width = Inf
)


readr::write_csv(
  categorical_tests,
  file.path(
    tables_dir,
    "clinical_categorical_tests_effect_sizes.csv"
  )
)

readr::write_csv(
  age_test,
  file.path(
    tables_dir,
    "clinical_age_test_effect_size.csv"
  )
)


# ============================================================
# PART B. SEX ASSOCIATIONS WITH THE FIRST THREE PCA SCORES
# ============================================================


# ------------------------------------------------------------
# 8. Recompute the baseline three-PC PCA
#
# The PCA is recomputed directly from X_snv so that this script is
# self-contained with respect to PCA scores. The specification is
# identical to the baseline analysis: no additional feature scaling
# and three retained PCs.
# ------------------------------------------------------------

pca_baseline <- FactoMineR::PCA(
  X_snv,
  scale.unit = FALSE,
  ncp = N_PC_BASELINE,
  graph = FALSE
)

scores_baseline <- as.matrix(
  pca_baseline$ind$coord[
    ,
    seq_len(N_PC_BASELINE),
    drop = FALSE
  ]
)

rownames(scores_baseline) <-
  sample_ids


scores_df <- tibble::tibble(
  participant_id = sample_ids,
  PC1 =
    scores_baseline[
      ,
      1
    ],
  PC2 =
    scores_baseline[
      ,
      2
    ],
  PC3 =
    scores_baseline[
      ,
      3
    ]
) |>
  dplyr::left_join(
    metadata |>
      dplyr::select(
        participant_id,
        sex,
        cluster
      ),
    by = "participant_id"
  )


# Effect-size sign convention used below:
#   rank-biserial correlation: positive = higher scores in F,
#                              negative = higher scores in M
#   Hedges' g:                F minus M
#   Hodges-Lehmann estimate:  F minus M
#
# Thus, negative values indicate higher scores in M than in F.
pc_names <- c(
  "PC1",
  "PC2",
  "PC3"
)


# Wilcoxon rank-sum tests provide the nonparametric group
# comparison. Rank-biserial correlation is the effect size most
# directly aligned with the rank-based test. Hedges' g and the
# Hodges-Lehmann location shift are also reported to aid interpretation.
pc_sex_results <- lapply(
  seq_along(pc_names),
  function(i) {

    pc_name <- pc_names[i]

    x <- scores_df[[pc_name]]

    sex <- scores_df$sex

    wilcox_res <- stats::wilcox.test(
      x ~ sex,
      conf.int = TRUE
    )

    x_f <- x[
      sex == "F"
    ]

    x_m <- x[
      sex == "M"
    ]

    r_rb <- rank_biserial_two_group(
      x_f,
      x_m
    )

    r_rb_ci <-
      bootstrap_two_group_effect(
        x = x,
        sex = sex,
        effect_fun =
          rank_biserial_two_group,
        B = N_BOOT_EFFECT,
        seed = SEED + 200L + i
      )

    g <- hedges_g_two_group(
      x_f,
      x_m
    )

    g_ci <-
      bootstrap_two_group_effect(
        x = x,
        sex = sex,
        effect_fun =
          hedges_g_two_group,
        B = N_BOOT_EFFECT,
        seed = SEED + 300L + i
      )

    tibble::tibble(
      PC = pc_name,
      n_F = length(x_f),
      n_M = length(x_m),
      mean_F = mean(x_f),
      mean_M = mean(x_m),
      median_F =
        stats::median(x_f),
      median_M =
        stats::median(x_m),
      wilcoxon_W =
        as.numeric(
          wilcox_res$statistic
        ),
      p_value =
        wilcox_res$p.value,
      hodges_lehmann_F_minus_M =
        as.numeric(
          wilcox_res$estimate
        ),
      hodges_lehmann_ci_lower =
        as.numeric(
          wilcox_res$conf.int[1]
        ),
      hodges_lehmann_ci_upper =
        as.numeric(
          wilcox_res$conf.int[2]
        ),
      rank_biserial_F_vs_M =
        r_rb,
      rank_biserial_ci_lower =
        r_rb_ci[
          "lower"
        ],
      rank_biserial_ci_upper =
        r_rb_ci[
          "upper"
        ],
      hedges_g_F_minus_M =
        g,
      hedges_g_ci_lower =
        g_ci[
          "lower"
        ],
      hedges_g_ci_upper =
        g_ci[
          "upper"
        ]
    )
  }
)

pc_sex_results <- dplyr::bind_rows(
  pc_sex_results
)


cat("\n========================================\n")
cat("PC1-PC3 SEX ASSOCIATIONS\n")
cat("========================================\n")

print(
  pc_sex_results,
  n = Inf,
  width = Inf
)


readr::write_csv(
  pc_sex_results,
  file.path(
    tables_dir,
    "clinical_PC1_PC3_sex_effect_sizes.csv"
  )
)


# ------------------------------------------------------------
# 9. PC1 visualization by sex
# ------------------------------------------------------------

p_pc1_sex <- ggplot2::ggplot(
  scores_df,
  ggplot2::aes(
    x = sex,
    y = PC1
  )
) +
  ggplot2::geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    height = 0,
    size = 2
  ) +
  ggplot2::labs(
    x = "Sex",
    y = "PC1 score",
    title = "Baseline PC1 scores by sex"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "clinical_PC1_by_sex.png"
  ),
  p_pc1_sex,
  width = 6,
  height = 5,
  dpi = 300
)


# ============================================================
# PART C. SEX-RESIDUALIZED SPECTRAL SENSITIVITY
# ============================================================
#
# For each of the 793 SNV-normalized spectral variables, the
# participant-level intensity is decomposed as:
#
#   intensity_j = intercept_j + beta_j(sex) + residual_j
#
# The residual matrix removes the sample-wide mean component
# associated with sex at each wavenumber while retaining the
# remaining participant-level spectral variation.
#
# Matrix least squares is used for efficiency; it is algebraically
# equivalent to fitting 793 separate ordinary least-squares models
# with the same design matrix (intercept + sex).
#
# This step is deliberately performed BEFORE recomputing PCA.
# Therefore both the PCA axes and the clustering solution are learned
# anew from the sex-residualized spectra rather than projecting
# residuals onto the original PCA space.
#
# PCA and clustering then use the same downstream specification as
# the baseline analysis: no additional variable scaling, three PCs,
# Euclidean distance, Ward.D2 linkage, and k = 3.
#
# ARI compares partitions independently of cluster-label numbering:
#   ARI = 1   identical partition
#   ARI ~ 0  agreement near chance
#   ARI < 0  less agreement than expected by chance
# ============================================================


# ------------------------------------------------------------
# 10. Remove the linear sex-associated component
# ------------------------------------------------------------

design_sex <- stats::model.matrix(
  ~ sex,
  data = metadata
)

coef_sex <- qr.solve(
  design_sex,
  X_snv
)

X_sex_residualized <-
  X_snv -
    design_sex %*%
    coef_sex

rownames(
  X_sex_residualized
) <- sample_ids

colnames(
  X_sex_residualized
) <- colnames(X_snv)


# Numerical checks: residuals should be orthogonal to the
# intercept and sex contrast up to floating-point tolerance.

residual_orthogonality <-
  crossprod(
    design_sex,
    X_sex_residualized
  )

max_abs_residual_orthogonality <-
  max(
    abs(
      residual_orthogonality
    )
  )


# ------------------------------------------------------------
# 11. Recompute PCA and k = 3 Ward.D2 clustering
# ------------------------------------------------------------

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

rownames(
  scores_sex_residualized
) <- sample_ids


cluster_sex_residualized_raw <-
  ward_clusters(
    scores_sex_residualized,
    k = K_BASELINE
  )

names(
  cluster_sex_residualized_raw
) <- sample_ids


# Cluster numbers are arbitrary. They are aligned to the baseline
# labels only to make participant-level tables readable and to count
# how many participants change an aligned cluster label. ARI itself
# is label-invariant.

cluster_sex_residualized <-
  align_cluster_labels(
    reference =
      cluster_baseline,
    candidate =
      cluster_sex_residualized_raw
  )


sex_residualized_ARI <-
  mclust::adjustedRandIndex(
    cluster_baseline,
    cluster_sex_residualized
  )


sex_residualized_sil <-
  cluster::silhouette(
    cluster_sex_residualized,
    stats::dist(
      scores_sex_residualized,
      method = "euclidean"
    )
  )

sex_residualized_average_silhouette <-
  mean(
    sex_residualized_sil[
      ,
      "sil_width"
    ]
  )


resid_eig <- as.data.frame(
  pca_sex_residualized$eig
)

resid_cumulative_variance_3PC <-
  resid_eig[
    N_PC_BASELINE,
    3
  ]


# ------------------------------------------------------------
# 12. Residual sex association after residualization
#
# This is a diagnostic check on the residualized partition. It asks
# whether cluster membership is still associated with sex after the
# linear sex-associated spectral component has been removed. It is
# not used to select the number of clusters.
# ------------------------------------------------------------

sex_resid_tab <- table(
  metadata$sex,
  factor(
    cluster_sex_residualized[
      metadata$participant_id
    ],
    levels = 1:K_BASELINE
  )
)

sex_resid_fisher <- stats::fisher.test(
  sex_resid_tab
)

sex_resid_v <- cramers_v(
  metadata$sex,
  factor(
    cluster_sex_residualized[
      metadata$participant_id
    ],
    levels = 1:K_BASELINE
  )
)

sex_resid_v_ci <- bootstrap_cramers_v(
  x = metadata$sex,
  y = factor(
    cluster_sex_residualized[
      metadata$participant_id
    ],
    levels = 1:K_BASELINE
  ),
  B = N_BOOT_EFFECT,
  seed = SEED + 500L
)


# ------------------------------------------------------------
# 13. Participant-level residualized assignments
# ------------------------------------------------------------

sex_residualized_assignments <- tibble::tibble(
  participant_id = sample_ids,
  sex =
    as.character(
      metadata$sex
    ),
  baseline_cluster =
    as.integer(
      cluster_baseline[
        sample_ids
      ]
    ),
  sex_residualized_cluster =
    as.integer(
      cluster_sex_residualized[
        sample_ids
      ]
    ),
  changed_cluster =
    baseline_cluster !=
      sex_residualized_cluster,
  residualized_PC1 =
    scores_sex_residualized[
      ,
      1
    ],
  residualized_PC2 =
    scores_sex_residualized[
      ,
      2
    ],
  residualized_PC3 =
    scores_sex_residualized[
      ,
      3
    ]
)


n_changed_clusters <- sum(
  sex_residualized_assignments$changed_cluster
)


sex_residualized_summary <- tibble::tibble(
  n_participants =
    nrow(X_snv),
  baseline_cluster_sizes =
    paste(
      as.integer(
        table(cluster_baseline)
      ),
      collapse = "/"
    ),
  residualized_cluster_sizes =
    paste(
      as.integer(
        table(
          factor(
            cluster_sex_residualized,
            levels = 1:K_BASELINE
          )
        )
      ),
      collapse = "/"
    ),
  ARI_vs_baseline =
    sex_residualized_ARI,
  n_changed_cluster_after_label_alignment =
    n_changed_clusters,
  residualized_average_silhouette =
    sex_residualized_average_silhouette,
  residualized_3PC_cumulative_variance_percent =
    resid_cumulative_variance_3PC,
  sex_cluster_fisher_p_after_residualization =
    sex_resid_fisher$p.value,
  sex_cluster_cramers_v_after_residualization =
    sex_resid_v,
  sex_cluster_cramers_v_ci_lower =
    sex_resid_v_ci[
      "lower"
    ],
  sex_cluster_cramers_v_ci_upper =
    sex_resid_v_ci[
      "upper"
    ],
  max_abs_design_residual_crossproduct =
    max_abs_residual_orthogonality
)


cat("\n========================================\n")
cat("SEX-RESIDUALIZED CLUSTERING SENSITIVITY\n")
cat("========================================\n")

print(
  sex_residualized_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  sex_residualized_assignments,
  file.path(
    tables_dir,
    "sex_residualized_cluster_assignments.csv"
  )
)

readr::write_csv(
  sex_residualized_summary,
  file.path(
    tables_dir,
    "sex_residualized_clustering_summary.csv"
  )
)

saveRDS(
  X_sex_residualized,
  file.path(
    processed_dir,
    "pediatric_X_snv_sex_residualized.rds"
  )
)

saveRDS(
  pca_sex_residualized,
  file.path(
    processed_dir,
    "pediatric_PCA_sex_residualized.rds"
  )
)

saveRDS(
  cluster_sex_residualized,
  file.path(
    processed_dir,
    "pediatric_clusters_sex_residualized.rds"
  )
)


# ------------------------------------------------------------
# 14. Baseline vs sex-residualized cluster cross-tabulation
# ------------------------------------------------------------

baseline_residualized_cross_tab <- as.data.frame.matrix(
  table(
    baseline_cluster =
      cluster_baseline[
        sample_ids
      ],
    sex_residualized_cluster =
      cluster_sex_residualized[
        sample_ids
      ]
  )
)

baseline_residualized_cross_tab <-
  tibble::rownames_to_column(
    baseline_residualized_cross_tab,
    var =
      "baseline_cluster"
  )


readr::write_csv(
  baseline_residualized_cross_tab,
  file.path(
    tables_dir,
    "baseline_vs_sex_residualized_clusters.csv"
  )
)


# ============================================================
# PART D. CONSOLIDATED REPRODUCIBILITY OUTPUTS
# ============================================================


# ------------------------------------------------------------
# 15. One compact metadata / sex-sensitivity summary
#
# This long-format table is the preferred numerical source for the
# Quarto report. Downstream text should read value_numeric directly
# rather than duplicate values manually.
# ------------------------------------------------------------

pc1_row <- pc_sex_results |>
  dplyr::filter(
    PC == "PC1"
  )

sex_baseline_row <- categorical_tests |>
  dplyr::filter(
    variable == "sex"
  )

age_row <- age_test


sex_sensitivity_summary <- dplyr::bind_rows(

  tibble::tibble(
    section =
      "Clinical metadata associations",
    metric = c(
      "Sex x cluster Fisher p",
      "Sex x cluster Cramer's V",
      "Sex x cluster Cramer's V CI lower",
      "Sex x cluster Cramer's V CI upper",
      "Age x cluster Kruskal-Wallis p",
      "Age epsilon-squared",
      "Endocrine x cluster Fisher p",
      "UCTD x cluster Fisher p",
      "Periodontal x cluster Fisher p"
    ),
    value_numeric = c(
      sex_baseline_row$p_value,
      sex_baseline_row$cramers_v,
      sex_baseline_row$cramers_v_ci_lower,
      sex_baseline_row$cramers_v_ci_upper,
      age_row$p_value,
      age_row$epsilon_squared,
      categorical_tests$p_value[
        categorical_tests$variable ==
          "endocrine_status"
      ],
      categorical_tests$p_value[
        categorical_tests$variable ==
          "uctd_status"
      ],
      categorical_tests$p_value[
        categorical_tests$variable ==
          "periodontal_status"
      ]
    ),
    value_text = format(
      c(
        sex_baseline_row$p_value,
        sex_baseline_row$cramers_v,
        sex_baseline_row$cramers_v_ci_lower,
        sex_baseline_row$cramers_v_ci_upper,
        age_row$p_value,
        age_row$epsilon_squared,
        categorical_tests$p_value[
          categorical_tests$variable ==
            "endocrine_status"
        ],
        categorical_tests$p_value[
          categorical_tests$variable ==
            "uctd_status"
        ],
        categorical_tests$p_value[
          categorical_tests$variable ==
            "periodontal_status"
        ]
      ),
      digits = 10
    )
  ),

  tibble::tibble(
    section =
      "PC1 sex association",
    metric = c(
      "Wilcoxon p",
      "Rank-biserial correlation F vs M",
      "Rank-biserial 95% CI lower",
      "Rank-biserial 95% CI upper",
      "Hedges g F minus M",
      "Hedges g 95% CI lower",
      "Hedges g 95% CI upper",
      "Hodges-Lehmann F minus M",
      "Hodges-Lehmann 95% CI lower",
      "Hodges-Lehmann 95% CI upper"
    ),
    value_numeric = c(
      pc1_row$p_value,
      pc1_row$rank_biserial_F_vs_M,
      pc1_row$rank_biserial_ci_lower,
      pc1_row$rank_biserial_ci_upper,
      pc1_row$hedges_g_F_minus_M,
      pc1_row$hedges_g_ci_lower,
      pc1_row$hedges_g_ci_upper,
      pc1_row$hodges_lehmann_F_minus_M,
      pc1_row$hodges_lehmann_ci_lower,
      pc1_row$hodges_lehmann_ci_upper
    ),
    value_text = format(
      c(
        pc1_row$p_value,
        pc1_row$rank_biserial_F_vs_M,
        pc1_row$rank_biserial_ci_lower,
        pc1_row$rank_biserial_ci_upper,
        pc1_row$hedges_g_F_minus_M,
        pc1_row$hedges_g_ci_lower,
        pc1_row$hedges_g_ci_upper,
        pc1_row$hodges_lehmann_F_minus_M,
        pc1_row$hodges_lehmann_ci_lower,
        pc1_row$hodges_lehmann_ci_upper
      ),
      digits = 10
    )
  ),

  tibble::tibble(
    section =
      "Sex-residualized sensitivity",
    metric = c(
      "ARI vs baseline",
      "Participants changing aligned cluster",
      "Residualized average silhouette",
      "Residualized 3-PC cumulative variance (%)",
      "Sex x residualized cluster Fisher p",
      "Sex x residualized cluster Cramer's V",
      "Sex x residualized cluster Cramer's V CI lower",
      "Sex x residualized cluster Cramer's V CI upper"
    ),
    value_numeric = c(
      sex_residualized_summary$ARI_vs_baseline,
      sex_residualized_summary$n_changed_cluster_after_label_alignment,
      sex_residualized_summary$residualized_average_silhouette,
      sex_residualized_summary$residualized_3PC_cumulative_variance_percent,
      sex_residualized_summary$sex_cluster_fisher_p_after_residualization,
      sex_residualized_summary$sex_cluster_cramers_v_after_residualization,
      sex_residualized_summary$sex_cluster_cramers_v_ci_lower,
      sex_residualized_summary$sex_cluster_cramers_v_ci_upper
    ),
    value_text = format(
      c(
        sex_residualized_summary$ARI_vs_baseline,
        sex_residualized_summary$n_changed_cluster_after_label_alignment,
        sex_residualized_summary$residualized_average_silhouette,
        sex_residualized_summary$residualized_3PC_cumulative_variance_percent,
        sex_residualized_summary$sex_cluster_fisher_p_after_residualization,
        sex_residualized_summary$sex_cluster_cramers_v_after_residualization,
        sex_residualized_summary$sex_cluster_cramers_v_ci_lower,
        sex_residualized_summary$sex_cluster_cramers_v_ci_upper
      ),
      digits = 10
    )
  )
)


readr::write_csv(
  sex_sensitivity_summary,
  file.path(
    tables_dir,
    "sex_sensitivity_summary.csv"
  )
)


# ------------------------------------------------------------
# 16. Save compact result bundle for Quarto / downstream reporting
#
# The bundle contains analysis-ready results rather than raw clinical
# notes. It allows the reproducible report to read numerical results
# programmatically without copying values by hand.
# ------------------------------------------------------------

sex_sensitivity_objects <- list(
  settings = list(
    seed = SEED,
    baseline_n_pc =
      N_PC_BASELINE,
    baseline_k =
      K_BASELINE,
    n_boot_effect =
      N_BOOT_EFFECT
  ),
  metadata =
    metadata_with_cluster,
  metadata_validation =
    metadata_validation,
  age_by_cluster =
    age_by_cluster,
  categorical_cluster_counts =
    categorical_cluster_counts,
  categorical_tests =
    categorical_tests,
  age_test =
    age_test,
  pc_sex_results =
    pc_sex_results,
  sex_residualized_assignments =
    sex_residualized_assignments,
  sex_residualized_summary =
    sex_residualized_summary,
  baseline_residualized_cross_tab =
    baseline_residualized_cross_tab,
  sex_sensitivity_summary =
    sex_sensitivity_summary
)


saveRDS(
  sex_sensitivity_objects,
  file.path(
    processed_dir,
    "sex_sensitivity_objects.rds"
  )
)


# ------------------------------------------------------------
# 17. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    tables_dir,
    "sessionInfo_sex_sensitivity.txt"
  )
)


# ------------------------------------------------------------
# 18. Final console report
# ------------------------------------------------------------

cat("\n\n")
cat("========================================\n")
cat("CLINICAL METADATA / SEX SENSITIVITY COMPLETE\n")
cat("========================================\n")

cat(
  "Sex x baseline cluster Fisher p: ",
  signif(
    categorical_tests$p_value[
      categorical_tests$variable == "sex"
    ],
    4
  ),
  "\n",
  sep = ""
)

cat(
  "PC1 sex Wilcoxon p: ",
  signif(
    pc_sex_results$p_value[
      pc_sex_results$PC == "PC1"
    ],
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Sex-residualized ARI vs baseline: ",
  round(
    sex_residualized_ARI,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "Residualized cluster sizes: ",
  sex_residualized_summary$residualized_cluster_sizes,
  "\n",
  sep = ""
)

cat(
  "Participants changing aligned cluster: ",
  n_changed_clusters,
  "\n",
  sep = ""
)

cat(
  "Sex x residualized cluster Fisher p: ",
  signif(
    sex_resid_fisher$p.value,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "\nKey reporting table:\n",
  normalizePath(
    file.path(
      tables_dir,
      "sex_sensitivity_summary.csv"
    ),
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "\nFull result bundle:\n",
  normalizePath(
    file.path(
      processed_dir,
      "sex_sensitivity_objects.rds"
    ),
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat("========================================\n")
