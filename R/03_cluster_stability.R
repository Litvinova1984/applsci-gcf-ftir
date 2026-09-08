# ============================================================
# 03_cluster_stability.R
#
# Pediatric GCF FTIR
# Robustness analysis of the baseline three-cluster solution
#
# Baseline analytical specification:
#   - input: row-wise SNV-normalized pediatric spectra
#   - PCA: no additional feature-wise scaling
#   - retained baseline dimensionality: 3 PCs
#   - clustering space: PC1-PC3 scores
#   - distance: Euclidean
#   - linkage: Ward.D2
#   - baseline number of clusters: k = 3
#
# Analyses:
#   A. Sensitivity to the number of retained PCs
#   B. Sensitivity to the number of clusters
#   C. Baseline cluster-wise silhouette diagnostics
#   C1. Additional internal characterization in the full
#       SNV spectral space (full-space silhouette, PERMANOVA,
#       and betadisper)
#   D. Leave-one-out stability
#   E. Repeated 80% subsampling and pairwise consensus
#   F. Participant-level bootstrap stability
#
# All numerical outputs used for reporting are written to
# tables/ and data/pediatric/processed/. In particular,
# baseline_internal_characterization_full_snv.csv keeps the
# full-spectral-space silhouette, PERMANOVA, and betadisper
# results separate from the primary 3-PC robustness analyses.
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
  "vegan",
  "tibble",
  "dplyr",
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
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


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

PC_VALUES <- 2:7
K_VALUES <- 2:6

N_SUBSAMPLE <- 1000L
SUBSAMPLE_FRACTION <- 0.80

N_BOOTSTRAP <- 1000L


# ------------------------------------------------------------
# 2. Load and validate baseline inputs
# ------------------------------------------------------------

X_snv <- as.matrix(
  readRDS(x_file)
)

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

sample_ids <- rownames(X_snv)

if (is.null(sample_ids)) {
  stop("X_snv must have participant IDs as row names.")
}


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

  cluster_baseline_raw <- cluster_baseline_raw[
    sample_ids
  ]

  if (is.factor(cluster_baseline_raw)) {
    cluster_baseline <- as.integer(
      as.character(cluster_baseline_raw)
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

baseline_sizes <- as.integer(
  table(cluster_baseline)
)

if (!identical(
  baseline_sizes,
  c(15L, 4L, 6L)
)) {
  stop(
    "Expected baseline cluster sizes 15 / 4 / 6; obtained ",
    paste(baseline_sizes, collapse = " / "),
    "."
  )
}


cat("\n========================================\n")
cat("BASELINE INPUT\n")
cat("========================================\n")

cat(
  "Spectra: ",
  nrow(X_snv),
  "\n",
  sep = ""
)

cat(
  "Variables: ",
  ncol(X_snv),
  "\n",
  sep = ""
)

cat(
  "Baseline PCs: ",
  N_PC_BASELINE,
  "\n",
  sep = ""
)

cat(
  "Baseline k: ",
  K_BASELINE,
  "\n",
  sep = ""
)

cat(
  "Baseline cluster sizes: ",
  paste(
    baseline_sizes,
    collapse = " / "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

ward_clusters <- function(scores, k) {

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


# Calinski-Harabasz index:
# larger values indicate better separated / more compact
# clustering for a fixed data representation.

calinski_harabasz <- function(X, cl) {

  X <- as.matrix(X)

  cl <- as.integer(
    as.factor(cl)
  )

  n <- nrow(X)
  k <- length(unique(cl))

  grand_centroid <- colMeans(X)

  within_ss <- 0
  between_ss <- 0

  for (g in sort(unique(cl))) {

    Xg <- X[
      cl == g,
      ,
      drop = FALSE
    ]

    centroid_g <- colMeans(Xg)

    within_ss <- within_ss +
      sum(
        rowSums(
          sweep(
            Xg,
            2,
            centroid_g,
            "-"
          )^2
        )
      )

    between_ss <- between_ss +
      nrow(Xg) *
      sum(
        (
          centroid_g -
            grand_centroid
        )^2
      )
  }

  (
    between_ss / (k - 1)
  ) /
    (
      within_ss / (n - k)
    )
}


# Davies-Bouldin index:
# smaller values indicate better clustering.

davies_bouldin <- function(X, cl) {

  X <- as.matrix(X)

  cl <- as.integer(
    as.factor(cl)
  )

  groups <- sort(
    unique(cl)
  )

  k <- length(groups)


  centroids <- do.call(
    rbind,
    lapply(
      groups,
      function(g) {

        colMeans(
          X[
            cl == g,
            ,
            drop = FALSE
          ]
        )
      }
    )
  )


  scatter <- vapply(
    groups,
    function(g) {

      Xg <- X[
        cl == g,
        ,
        drop = FALSE
      ]

      cg <- colMeans(Xg)

      mean(
        sqrt(
          rowSums(
            sweep(
              Xg,
              2,
              cg,
              "-"
            )^2
          )
        )
      )
    },
    numeric(1)
  )


  centroid_dist <- as.matrix(
    stats::dist(
      centroids,
      method = "euclidean"
    )
  )


  R <- matrix(
    NA_real_,
    nrow = k,
    ncol = k
  )


  for (i in seq_len(k)) {

    for (j in seq_len(k)) {

      if (i != j) {

        R[i, j] <- (
          scatter[i] +
            scatter[j]
        ) /
          centroid_dist[i, j]
      }
    }
  }


  mean(
    apply(
      R,
      1,
      max,
      na.rm = TRUE
    )
  )
}


# Cluster-specific bootstrap stability:
# for each baseline cluster, identify the candidate cluster
# with maximum Jaccard membership overlap. Cluster labels are
# arbitrary and are therefore not compared directly.

max_cluster_jaccard <- function(
    reference,
    candidate,
    reference_cluster
) {

  reference_members <- names(reference)[
    reference == reference_cluster
  ]

  candidate_labels <- sort(
    unique(candidate)
  )

  jaccards <- vapply(
    candidate_labels,
    function(g) {

      candidate_members <- names(candidate)[
        candidate == g
      ]

      intersection_n <- length(
        intersect(
          reference_members,
          candidate_members
        )
      )

      union_n <- length(
        union(
          reference_members,
          candidate_members
        )
      )

      intersection_n / union_n
    },
    numeric(1)
  )

  max(jaccards)
}


# ------------------------------------------------------------
# 4. Recompute baseline PCA and verify baseline partition
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

rownames(scores_baseline) <- sample_ids


cluster_recomputed <- ward_clusters(
  scores_baseline,
  k = K_BASELINE
)

baseline_reproduction_ARI <-
  mclust::adjustedRandIndex(
    cluster_baseline,
    cluster_recomputed
  )

if (
  abs(
    baseline_reproduction_ARI - 1
  ) > 1e-12
) {
  stop(
    "The baseline three-PC Ward.D2 solution was not ",
    "reproduced exactly."
  )
}


baseline_eig <- as.data.frame(
  pca_baseline$eig
)

baseline_cumulative_variance <-
  baseline_eig[
    N_PC_BASELINE,
    3
  ]


# ============================================================
# A. PCA-DIMENSION SENSITIVITY
# ============================================================
#
# The baseline uses three PCs. Alternative dimensionalities
# from two to seven PCs are evaluated only as a downstream
# robustness analysis. They do not redefine the baseline
# PCA-retention decision.
#
# For each dimensionality:
#   - PCA is recalculated from X_snv;
#   - Ward.D2 clustering is repeated at k = 3;
#   - agreement with the three-PC baseline is quantified
#     using the Adjusted Rand Index (ARI).
# ============================================================

pc_results <- lapply(
  PC_VALUES,
  function(n_pc) {

    pca_i <- FactoMineR::PCA(
      X_snv,
      scale.unit = FALSE,
      ncp = n_pc,
      graph = FALSE
    )

    scores_i <- as.matrix(
      pca_i$ind$coord[
        ,
        seq_len(n_pc),
        drop = FALSE
      ]
    )

    cl_i <- ward_clusters(
      scores_i,
      k = K_BASELINE
    )

    ari_i <- mclust::adjustedRandIndex(
      cluster_baseline,
      cl_i
    )

    eig_i <- as.data.frame(
      pca_i$eig
    )

    tibble::tibble(
      n_PC = n_pc,
      cumulative_variance_percent =
        eig_i[n_pc, 3],
      ARI_vs_baseline_3PC = ari_i,
      identical_to_baseline =
        abs(ari_i - 1) < 1e-12,
      cluster_sizes =
        paste(
          as.integer(
            table(cl_i)
          ),
          collapse = "/"
        )
    )
  }
)

pc_results <- dplyr::bind_rows(
  pc_results
)


cat("\n========================================\n")
cat("PCA-DIMENSION SENSITIVITY: 2-7 PCs\n")
cat("========================================\n")

print(
  pc_results,
  n = Inf,
  width = Inf
)


readr::write_csv(
  pc_results,
  file.path(
    tables_dir,
    "cluster_PC_count_sensitivity.csv"
  )
)


p_pc_sensitivity <- ggplot2::ggplot(
  pc_results,
  ggplot2::aes(
    x = n_PC,
    y = ARI_vs_baseline_3PC
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(
    size = 2.5
  ) +
  ggplot2::geom_vline(
    xintercept = N_PC_BASELINE,
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dotted"
  ) +
  ggplot2::scale_x_continuous(
    breaks = PC_VALUES
  ) +
  ggplot2::labs(
    x = "Number of retained principal components",
    y = "Adjusted Rand Index vs 3-PC baseline",
    title = "Sensitivity of clustering to PCA dimensionality"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "cluster_PC_count_sensitivity.png"
  ),
  p_pc_sensitivity,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# B. CLUSTER-NUMBER SENSITIVITY
# ============================================================
#
# Alternative k = 2:6 solutions are evaluated in the fixed
# baseline three-PC score space. Three complementary internal
# validity measures are reported:
#
#   - average silhouette: larger is better;
#   - Calinski-Harabasz: larger is better;
#   - Davies-Bouldin: smaller is better.
#
# Cluster sizes are retained because apparent improvements
# at larger k can coincide with fragmentation into very small
# clusters in a small pilot sample.
# ============================================================

k_results <- lapply(
  K_VALUES,
  function(k) {

    cl <- ward_clusters(
      scores_baseline,
      k = k
    )

    d <- stats::dist(
      scores_baseline,
      method = "euclidean"
    )

    sil <- cluster::silhouette(
      cl,
      d
    )

    size_vector <- as.integer(
      table(cl)
    )

    tibble::tibble(
      k = k,
      average_silhouette =
        mean(
          sil[, "sil_width"]
        ),
      calinski_harabasz =
        calinski_harabasz(
          scores_baseline,
          cl
        ),
      davies_bouldin =
        davies_bouldin(
          scores_baseline,
          cl
        ),
      smallest_cluster =
        min(size_vector),
      largest_cluster =
        max(size_vector),
      cluster_sizes =
        paste(
          size_vector,
          collapse = "/"
        )
    )
  }
)

k_results <- dplyr::bind_rows(
  k_results
)


cat("\n========================================\n")
cat("CLUSTER-NUMBER SENSITIVITY: k = 2-6\n")
cat("========================================\n")

print(
  k_results,
  n = Inf,
  width = Inf
)


readr::write_csv(
  k_results,
  file.path(
    tables_dir,
    "cluster_number_sensitivity_k2_k6.csv"
  )
)


# Focused comparison of the parsimonious k = 2 and k = 3
# solutions. This table is convenient for manuscript/rebuttal
# reporting; the full k = 2:6 table remains the primary record.

k2_row <- k_results |>
  dplyr::filter(
    k == 2
  )

k3_row <- k_results |>
  dplyr::filter(
    k == K_BASELINE
  )

k3_vs_k2 <- tibble::tibble(
  metric = c(
    "Average silhouette",
    "Calinski-Harabasz",
    "Davies-Bouldin"
  ),
  k2 = c(
    k2_row$average_silhouette,
    k2_row$calinski_harabasz,
    k2_row$davies_bouldin
  ),
  k3 = c(
    k3_row$average_silhouette,
    k3_row$calinski_harabasz,
    k3_row$davies_bouldin
  ),
  k3_minus_k2 = c(
    k3_row$average_silhouette -
      k2_row$average_silhouette,
    k3_row$calinski_harabasz -
      k2_row$calinski_harabasz,
    k3_row$davies_bouldin -
      k2_row$davies_bouldin
  )
)

readr::write_csv(
  k3_vs_k2,
  file.path(
    tables_dir,
    "cluster_number_k3_vs_k2.csv"
  )
)


p_silhouette <- ggplot2::ggplot(
  k_results,
  ggplot2::aes(
    x = k,
    y = average_silhouette
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_vline(
    xintercept = K_BASELINE,
    linetype = "dashed"
  ) +
  ggplot2::scale_x_continuous(
    breaks = K_VALUES
  ) +
  ggplot2::labs(
    x = "Number of clusters (k)",
    y = "Average silhouette width",
    title = "Cluster-number sensitivity"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "cluster_number_silhouette.png"
  ),
  p_silhouette,
  width = 7,
  height = 5,
  dpi = 300
)


p_ch <- ggplot2::ggplot(
  k_results,
  ggplot2::aes(
    x = k,
    y = calinski_harabasz
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_vline(
    xintercept = K_BASELINE,
    linetype = "dashed"
  ) +
  ggplot2::scale_x_continuous(
    breaks = K_VALUES
  ) +
  ggplot2::labs(
    x = "Number of clusters (k)",
    y = "Calinski-Harabasz index",
    title = "Cluster-number sensitivity"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "cluster_number_calinski_harabasz.png"
  ),
  p_ch,
  width = 7,
  height = 5,
  dpi = 300
)


p_db <- ggplot2::ggplot(
  k_results,
  ggplot2::aes(
    x = k,
    y = davies_bouldin
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_vline(
    xintercept = K_BASELINE,
    linetype = "dashed"
  ) +
  ggplot2::scale_x_continuous(
    breaks = K_VALUES
  ) +
  ggplot2::labs(
    x = "Number of clusters (k)",
    y = "Davies-Bouldin index",
    title = "Cluster-number sensitivity"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "cluster_number_davies_bouldin.png"
  ),
  p_db,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# C. BASELINE CLUSTER-WISE SILHOUETTE
# ============================================================
#
# The overall silhouette is complemented by cluster-specific
# summaries and participant-level widths. Negative widths
# indicate observations closer to another cluster than to
# their assigned cluster.
# ============================================================

baseline_dist <- stats::dist(
  scores_baseline,
  method = "euclidean"
)

baseline_sil <- cluster::silhouette(
  cluster_baseline,
  baseline_dist
)

baseline_sil_df <- tibble::tibble(
  sample_id = rownames(baseline_sil),
  cluster = as.integer(
    baseline_sil[, "cluster"]
  ),
  silhouette_width = as.numeric(
    baseline_sil[, "sil_width"]
  )
)

baseline_average_silhouette <- mean(
  baseline_sil_df$silhouette_width
)


cluster_silhouette_summary <- baseline_sil_df |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_silhouette =
      mean(silhouette_width),
    median_silhouette =
      stats::median(silhouette_width),
    min_silhouette =
      min(silhouette_width),
    max_silhouette =
      max(silhouette_width),
    n_negative =
      sum(silhouette_width < 0),
    proportion_negative =
      mean(silhouette_width < 0),
    .groups = "drop"
  )


cat("\n========================================\n")
cat("BASELINE CLUSTER-WISE SILHOUETTE\n")
cat("========================================\n")

cat(
  "Overall average silhouette: ",
  baseline_average_silhouette,
  "\n",
  sep = ""
)

print(
  cluster_silhouette_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  baseline_sil_df,
  file.path(
    tables_dir,
    "baseline_silhouette_by_participant.csv"
  )
)

readr::write_csv(
  cluster_silhouette_summary,
  file.path(
    tables_dir,
    "baseline_silhouette_by_cluster.csv"
  )
)


# ============================================================
# C1. ADDITIONAL INTERNAL CHARACTERIZATION IN FULL SNV SPACE
# ============================================================
#
# The baseline clustering itself is defined in the retained
# three-PC score space, so the 3-PC silhouette above is the
# geometry directly matched to the clustering procedure.
#
# For continuity with the original spectral-space description,
# this section also evaluates the same fixed cluster labels in
# the complete 793-variable SNV-normalized spectral matrix.
#
# These quantities answer a different question from the
# resampling and k-sensitivity analyses above:
#
#   - full-space silhouette describes separation when distances
#     are calculated from all retained SNV spectral variables;
#   - PERMANOVA quantifies how much of the total spectral
#     variation is associated with the already-defined clusters;
#   - betadisper checks whether between-cluster separation could
#     be accompanied by unequal within-cluster dispersion.
#
# Because cluster labels were derived from the same spectra,
# PERMANOVA and betadisper are treated as internal
# characterization, not as independent validation of the
# existence or biological meaning of the clusters.
#
# The full-space and 3-PC silhouette values are not expected to
# be identical because they use different distance geometries.
# ============================================================

full_snv_dist <- stats::dist(
  X_snv,
  method = "euclidean"
)

full_snv_sil <- cluster::silhouette(
  cluster_baseline,
  full_snv_dist
)

full_snv_sil_df <- tibble::tibble(
  sample_id = rownames(full_snv_sil),
  cluster = as.integer(
    full_snv_sil[, "cluster"]
  ),
  silhouette_width = as.numeric(
    full_snv_sil[, "sil_width"]
  )
)

full_snv_average_silhouette <- mean(
  full_snv_sil_df$silhouette_width
)


# PERMANOVA in the complete SNV spectral space.
# The pseudo-F statistic and R2 do not depend on the random
# permutation seed; the permutation p-value does.

cluster_factor <- factor(
  cluster_baseline,
  levels = sort(
    unique(cluster_baseline)
  )
)

set.seed(SEED)

permanova_full_snv <- vegan::adonis2(
  X_snv ~ cluster_factor,
  method = "euclidean",
  permutations = 9999,
  by = "terms"
)

# With by = "terms", the predictor is returned explicitly as
# the "cluster_factor" row. This makes extraction unambiguous
# and avoids relying on the default omnibus row name ("Model").
permanova_cluster_row <-
  permanova_full_snv[
    "cluster_factor",
    ,
    drop = FALSE
  ]

if (
  nrow(permanova_cluster_row) != 1L ||
  anyNA(
    c(
      permanova_cluster_row$R2,
      permanova_cluster_row$F,
      permanova_cluster_row$`Pr(>F)`
    )
  )
) {
  stop(
    "Could not extract the cluster term from the PERMANOVA table."
  )
}

permanova_R2 <- as.numeric(
  permanova_cluster_row$R2
)

permanova_F <- as.numeric(
  permanova_cluster_row$F
)

permanova_p <- as.numeric(
  permanova_cluster_row$`Pr(>F)`
)


# Homogeneity of multivariate dispersion in the same complete
# SNV spectral space. Both the classical ANOVA result and the
# permutation result are retained; the permutation p-value is
# the preferred quantity for reporting alongside PERMANOVA.

betadisper_full_snv <- vegan::betadisper(
  full_snv_dist,
  group = cluster_factor
)

betadisper_anova <- stats::anova(
  betadisper_full_snv
)

betadisper_anova_F <- as.numeric(
  betadisper_anova[
    1,
    "F value"
  ]
)

betadisper_anova_p <- as.numeric(
  betadisper_anova[
    1,
    "Pr(>F)"
  ]
)

set.seed(SEED)

betadisper_perm <- vegan::permutest(
  betadisper_full_snv,
  permutations = 9999
)

betadisper_perm_F <- as.numeric(
  betadisper_perm$tab[
    1,
    "F"
  ]
)

betadisper_perm_p <- as.numeric(
  betadisper_perm$tab[
    1,
    "Pr(>F)"
  ]
)


baseline_internal_characterization <-
  tibble::tibble(
    metric = c(
      "Average silhouette in PC1-PC3 clustering space",
      "Average silhouette in full SNV spectral space",
      "PERMANOVA R2 in full SNV spectral space",
      "PERMANOVA pseudo-F in full SNV spectral space",
      "PERMANOVA permutation p in full SNV spectral space",
      "Betadisper ANOVA F in full SNV spectral space",
      "Betadisper ANOVA p in full SNV spectral space",
      "Betadisper permutation F in full SNV spectral space",
      "Betadisper permutation p in full SNV spectral space"
    ),
    value = c(
      baseline_average_silhouette,
      full_snv_average_silhouette,
      permanova_R2,
      permanova_F,
      permanova_p,
      betadisper_anova_F,
      betadisper_anova_p,
      betadisper_perm_F,
      betadisper_perm_p
    )
  )


cat("\n========================================\n")
cat("ADDITIONAL FULL-SNV INTERNAL CHARACTERIZATION\n")
cat("========================================\n")

print(
  baseline_internal_characterization,
  n = Inf,
  width = Inf
)


readr::write_csv(
  full_snv_sil_df,
  file.path(
    tables_dir,
    "baseline_silhouette_full_snv_by_participant.csv"
  )
)

readr::write_csv(
  baseline_internal_characterization,
  file.path(
    tables_dir,
    "baseline_internal_characterization_full_snv.csv"
  )
)


# ============================================================
# D. LEAVE-ONE-OUT STABILITY
# ============================================================
#
# Each participant is removed in turn. PCA is recomputed from
# the remaining 24 participants, three PCs are retained, and
# Ward.D2 clustering is repeated at k = 3. The resulting
# partition is compared with the baseline partition restricted
# to the same 24 participants.
# ============================================================

loo_results <- lapply(
  sample_ids,
  function(removed_id) {

    keep_ids <- setdiff(
      sample_ids,
      removed_id
    )

    X_loo <- X_snv[
      keep_ids,
      ,
      drop = FALSE
    ]

    pca_loo <- FactoMineR::PCA(
      X_loo,
      scale.unit = FALSE,
      ncp = N_PC_BASELINE,
      graph = FALSE
    )

    scores_loo <- as.matrix(
      pca_loo$ind$coord[
        ,
        seq_len(N_PC_BASELINE),
        drop = FALSE
      ]
    )

    cl_loo <- ward_clusters(
      scores_loo,
      k = K_BASELINE
    )

    baseline_loo <- cluster_baseline[
      keep_ids
    ]

    ari <- mclust::adjustedRandIndex(
      baseline_loo,
      cl_loo
    )

    tibble::tibble(
      removed_sample = removed_id,
      ARI = ari,
      cluster_sizes =
        paste(
          as.integer(
            table(cl_loo)
          ),
          collapse = "/"
        )
    )
  }
)

loo_results <- dplyr::bind_rows(
  loo_results
)


loo_summary <- tibble::tibble(
  n_leave_one_out_runs =
    nrow(loo_results),
  mean_ARI =
    mean(loo_results$ARI),
  median_ARI =
    stats::median(loo_results$ARI),
  min_ARI =
    min(loo_results$ARI),
  max_ARI =
    max(loo_results$ARI),
  n_ARI_1 =
    sum(
      abs(
        loo_results$ARI - 1
      ) < 1e-12
    ),
  proportion_ARI_1 =
    mean(
      abs(
        loo_results$ARI - 1
      ) < 1e-12
    ),
  n_ARI_ge_0_8 =
    sum(
      loo_results$ARI >= 0.8
    ),
  proportion_ARI_ge_0_8 =
    mean(
      loo_results$ARI >= 0.8
    )
)


loo_influence <- loo_results |>
  dplyr::mutate(
    baseline_cluster =
      cluster_baseline[
        removed_sample
      ]
  ) |>
  dplyr::arrange(
    ARI
  )


cat("\n========================================\n")
cat("LEAVE-ONE-OUT STABILITY\n")
cat("========================================\n")

print(
  loo_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  loo_results,
  file.path(
    tables_dir,
    "cluster_LOO_results.csv"
  )
)

readr::write_csv(
  loo_summary,
  file.path(
    tables_dir,
    "cluster_LOO_summary.csv"
  )
)

readr::write_csv(
  loo_influence,
  file.path(
    tables_dir,
    "cluster_LOO_influence.csv"
  )
)


# ============================================================
# E. REPEATED 80% SUBSAMPLING AND CONSENSUS
# ============================================================
#
# In each of 1000 iterations:
#   - 80% of participants are sampled without replacement;
#   - PCA is recomputed;
#   - the baseline three PCs are retained;
#   - Ward.D2 clustering is repeated at k = 3;
#   - ARI is calculated against the baseline partition
#     restricted to the same participants.
#
# Pairwise consensus is the proportion of iterations in which
# two participants are assigned to the same cluster among
# iterations in which both participants are sampled.
# ============================================================

SUBSAMPLE_N <- floor(
  nrow(X_snv) *
    SUBSAMPLE_FRACTION
)

set.seed(SEED)


subsample_ari <- numeric(
  N_SUBSAMPLE
)

co_sampled <- matrix(
  0L,
  nrow = nrow(X_snv),
  ncol = nrow(X_snv),
  dimnames = list(
    sample_ids,
    sample_ids
  )
)

co_clustered <- matrix(
  0L,
  nrow = nrow(X_snv),
  ncol = nrow(X_snv),
  dimnames = list(
    sample_ids,
    sample_ids
  )
)


for (b in seq_len(N_SUBSAMPLE)) {

  selected_ids <- sample(
    sample_ids,
    size = SUBSAMPLE_N,
    replace = FALSE
  )

  X_b <- X_snv[
    selected_ids,
    ,
    drop = FALSE
  ]

  pca_b <- FactoMineR::PCA(
    X_b,
    scale.unit = FALSE,
    ncp = N_PC_BASELINE,
    graph = FALSE
  )

  scores_b <- as.matrix(
    pca_b$ind$coord[
      ,
      seq_len(N_PC_BASELINE),
      drop = FALSE
    ]
  )

  cl_b <- ward_clusters(
    scores_b,
    k = K_BASELINE
  )

  baseline_b <- cluster_baseline[
    selected_ids
  ]

  subsample_ari[b] <-
    mclust::adjustedRandIndex(
      baseline_b,
      cl_b
    )


  co_sampled[
    selected_ids,
    selected_ids
  ] <- co_sampled[
    selected_ids,
    selected_ids
  ] + 1L


  same_cluster_b <- outer(
    cl_b,
    cl_b,
    FUN = "=="
  )

  storage.mode(
    same_cluster_b
  ) <- "integer"


  co_clustered[
    selected_ids,
    selected_ids
  ] <- co_clustered[
    selected_ids,
    selected_ids
  ] +
    same_cluster_b
}


subsampling_results <- tibble::tibble(
  iteration = seq_len(
    N_SUBSAMPLE
  ),
  ARI = subsample_ari
)


subsampling_summary <- tibble::tibble(
  iterations = N_SUBSAMPLE,
  sampling_fraction =
    SUBSAMPLE_FRACTION,
  participants_per_iteration =
    SUBSAMPLE_N,
  mean_ARI =
    mean(subsample_ari),
  median_ARI =
    stats::median(subsample_ari),
  sd_ARI =
    stats::sd(subsample_ari),
  min_ARI =
    min(subsample_ari),
  q025_ARI =
    as.numeric(
      stats::quantile(
        subsample_ari,
        0.025
      )
    ),
  q25_ARI =
    as.numeric(
      stats::quantile(
        subsample_ari,
        0.25
      )
    ),
  q75_ARI =
    as.numeric(
      stats::quantile(
        subsample_ari,
        0.75
      )
    ),
  q975_ARI =
    as.numeric(
      stats::quantile(
        subsample_ari,
        0.975
      )
    ),
  max_ARI =
    max(subsample_ari),
  n_ARI_1 =
    sum(
      abs(
        subsample_ari - 1
      ) < 1e-12
    ),
  proportion_ARI_1 =
    mean(
      abs(
        subsample_ari - 1
      ) < 1e-12
    ),
  n_ARI_ge_0_8 =
    sum(
      subsample_ari >= 0.8
    ),
  proportion_ARI_ge_0_8 =
    mean(
      subsample_ari >= 0.8
    )
)


cat("\n========================================\n")
cat("REPEATED 80% SUBSAMPLING\n")
cat("========================================\n")

print(
  subsampling_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  subsampling_results,
  file.path(
    tables_dir,
    "cluster_subsampling_ARI_1000.csv"
  )
)

readr::write_csv(
  subsampling_summary,
  file.path(
    tables_dir,
    "cluster_subsampling_summary.csv"
  )
)


# ------------------------------------------------------------
# E1. Pairwise consensus matrix
# ------------------------------------------------------------

consensus_matrix <- matrix(
  NA_real_,
  nrow = nrow(co_sampled),
  ncol = ncol(co_sampled),
  dimnames = dimnames(co_sampled)
)

valid_pairs <- co_sampled > 0

consensus_matrix[
  valid_pairs
] <- co_clustered[
  valid_pairs
] /
  co_sampled[
    valid_pairs
  ]

diag(
  consensus_matrix
) <- 1


saveRDS(
  consensus_matrix,
  file.path(
    processed_dir,
    "cluster_consensus_matrix_80pct.rds"
  )
)

readr::write_csv(
  data.frame(
    sample_id =
      rownames(
        consensus_matrix
      ),
    consensus_matrix,
    check.names = FALSE
  ),
  file.path(
    tables_dir,
    "cluster_consensus_matrix_80pct.csv"
  )
)


pair_indices <- which(
  upper.tri(consensus_matrix),
  arr.ind = TRUE
)

pair_consensus <- tibble::tibble(
  sample_1 =
    rownames(consensus_matrix)[
      pair_indices[, 1]
    ],
  sample_2 =
    colnames(consensus_matrix)[
      pair_indices[, 2]
    ],
  baseline_same_cluster =
    cluster_baseline[
      sample_1
    ] ==
      cluster_baseline[
        sample_2
      ],
  consensus =
    consensus_matrix[
      pair_indices
    ]
)


consensus_summary <- pair_consensus |>
  dplyr::group_by(
    baseline_same_cluster
  ) |>
  dplyr::summarise(
    n_pairs =
      dplyr::n(),
    mean_consensus =
      mean(
        consensus,
        na.rm = TRUE
      ),
    median_consensus =
      stats::median(
        consensus,
        na.rm = TRUE
      ),
    min_consensus =
      min(
        consensus,
        na.rm = TRUE
      ),
    max_consensus =
      max(
        consensus,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


cat("\n========================================\n")
cat("PAIRWISE CONSENSUS SUMMARY\n")
cat("========================================\n")

print(
  consensus_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  pair_consensus,
  file.path(
    tables_dir,
    "cluster_consensus_pairwise.csv"
  )
)

readr::write_csv(
  consensus_summary,
  file.path(
    tables_dir,
    "cluster_consensus_summary.csv"
  )
)


# ------------------------------------------------------------
# E2. Consensus heatmap
# ------------------------------------------------------------

baseline_order <- order(
  cluster_baseline,
  names(cluster_baseline)
)

ordered_ids <- names(
  cluster_baseline
)[
  baseline_order
]

consensus_ordered <- consensus_matrix[
  ordered_ids,
  ordered_ids,
  drop = FALSE
]

heatmap_df <- as.data.frame(
  as.table(
    consensus_ordered
  )
)

colnames(
  heatmap_df
) <- c(
  "sample_1",
  "sample_2",
  "consensus"
)


p_consensus <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(
    x = sample_1,
    y = sample_2,
    fill = consensus
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    na.value = "grey90"
  ) +
  ggplot2::labs(
    x = "Participant",
    y = "Participant",
    fill = "Consensus",
    title = "Consensus clustering under repeated 80% subsampling"
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
    panel.grid =
      ggplot2::element_blank()
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "cluster_consensus_heatmap_80pct.png"
  ),
  p_consensus,
  width = 8,
  height = 7,
  dpi = 300
)


p_subsample_ari <- ggplot2::ggplot(
  subsampling_results,
  ggplot2::aes(
    x = ARI
  )
) +
  ggplot2::geom_histogram(
    bins = 30
  ) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  ggplot2::labs(
    x = "Adjusted Rand Index vs baseline",
    y = "Number of subsamples",
    title = "Cluster stability under repeated 80% subsampling"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "cluster_subsampling_ARI_distribution.png"
  ),
  p_subsample_ari,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# F. PARTICIPANT-LEVEL BOOTSTRAP STABILITY
# ============================================================
#
# Non-parametric participant-level bootstrap.
#
# In each of 1000 iterations:
#   1. sample n participants with replacement;
#   2. recompute PCA on the bootstrap sample;
#   3. retain the baseline three PCs;
#   4. perform Ward.D2 clustering at k = 3;
#   5. project all original participants into that bootstrap
#      PCA space;
#   6. assign each original participant to the nearest
#      bootstrap-cluster centroid;
#   7. compare the induced full-cohort partition with the
#      baseline partition using ARI.
#
# Cluster-specific Jaccard similarity is also reported.
# ============================================================

set.seed(SEED)


bootstrap_results <- vector(
  "list",
  N_BOOTSTRAP
)

bootstrap_jaccard <- vector(
  "list",
  N_BOOTSTRAP
)


for (b in seq_len(N_BOOTSTRAP)) {

  boot_ids <- sample(
    sample_ids,
    size = length(sample_ids),
    replace = TRUE
  )


  X_boot <- X_snv[
    boot_ids,
    ,
    drop = FALSE
  ]


  # The same participant may occur more than once in a
  # bootstrap sample. Unique row names are required by PCA.

  rownames(X_boot) <- paste0(
    "boot_",
    b,
    "_",
    seq_len(nrow(X_boot))
  )


  pca_boot <- FactoMineR::PCA(
    X_boot,
    scale.unit = FALSE,
    ncp = N_PC_BASELINE,
    graph = FALSE
  )


  scores_boot <- as.matrix(
    pca_boot$ind$coord[
      ,
      seq_len(N_PC_BASELINE),
      drop = FALSE
    ]
  )


  cl_boot <- ward_clusters(
    scores_boot,
    k = K_BASELINE
  )


  boot_cluster_labels <- sort(
    unique(cl_boot)
  )


  centroids_boot <- do.call(
    rbind,
    lapply(
      boot_cluster_labels,
      function(g) {

        colMeans(
          scores_boot[
            cl_boot == g,
            ,
            drop = FALSE
          ]
        )
      }
    )
  )


  rownames(centroids_boot) <-
    as.character(
      boot_cluster_labels
    )


  # Project all original participants into the PCA model
  # fitted to the current bootstrap sample.

  projection <- predict(
    pca_boot,
    newdata = as.data.frame(
      X_snv
    )
  )


  scores_all <- as.matrix(
    projection$coord[
      ,
      seq_len(N_PC_BASELINE),
      drop = FALSE
    ]
  )

  rownames(scores_all) <- sample_ids


  distance_to_centroids <- sapply(
    seq_len(nrow(centroids_boot)),
    function(j) {

      rowSums(
        sweep(
          scores_all,
          2,
          centroids_boot[j, ],
          "-"
        )^2
      )
    }
  )


  if (is.null(dim(distance_to_centroids))) {
    distance_to_centroids <- matrix(
      distance_to_centroids,
      ncol = 1
    )
  }


  nearest_centroid <- max.col(
    -distance_to_centroids,
    ties.method = "first"
  )


  induced_clusters <-
    boot_cluster_labels[
      nearest_centroid
    ]

  names(induced_clusters) <- sample_ids


  ari_b <- mclust::adjustedRandIndex(
    cluster_baseline,
    induced_clusters
  )


  bootstrap_results[[b]] <- tibble::tibble(
    iteration = b,
    ARI = ari_b,
    unique_participants =
      length(
        unique(boot_ids)
      ),
    induced_cluster_sizes =
      paste(
        as.integer(
          table(induced_clusters)
        ),
        collapse = "/"
      )
  )


  bootstrap_jaccard[[b]] <- tibble::tibble(
    iteration = b,
    baseline_cluster =
      sort(
        unique(cluster_baseline)
      ),
    Jaccard = vapply(
      sort(
        unique(cluster_baseline)
      ),
      function(g) {

        max_cluster_jaccard(
          reference =
            cluster_baseline,
          candidate =
            induced_clusters,
          reference_cluster =
            g
        )
      },
      numeric(1)
    )
  )
}


bootstrap_results <- dplyr::bind_rows(
  bootstrap_results
)

bootstrap_jaccard <- dplyr::bind_rows(
  bootstrap_jaccard
)


bootstrap_summary <- tibble::tibble(
  iterations = N_BOOTSTRAP,
  mean_ARI =
    mean(
      bootstrap_results$ARI
    ),
  median_ARI =
    stats::median(
      bootstrap_results$ARI
    ),
  sd_ARI =
    stats::sd(
      bootstrap_results$ARI
    ),
  min_ARI =
    min(
      bootstrap_results$ARI
    ),
  q025_ARI =
    as.numeric(
      stats::quantile(
        bootstrap_results$ARI,
        0.025
      )
    ),
  q25_ARI =
    as.numeric(
      stats::quantile(
        bootstrap_results$ARI,
        0.25
      )
    ),
  q75_ARI =
    as.numeric(
      stats::quantile(
        bootstrap_results$ARI,
        0.75
      )
    ),
  q975_ARI =
    as.numeric(
      stats::quantile(
        bootstrap_results$ARI,
        0.975
      )
    ),
  max_ARI =
    max(
      bootstrap_results$ARI
    ),
  n_ARI_1 =
    sum(
      abs(
        bootstrap_results$ARI - 1
      ) < 1e-12
    ),
  proportion_ARI_1 =
    mean(
      abs(
        bootstrap_results$ARI - 1
      ) < 1e-12
    ),
  n_ARI_ge_0_8 =
    sum(
      bootstrap_results$ARI >= 0.8
    ),
  proportion_ARI_ge_0_8 =
    mean(
      bootstrap_results$ARI >= 0.8
    ),
  mean_unique_participants =
    mean(
      bootstrap_results$unique_participants
    )
)


bootstrap_jaccard_summary <-
  bootstrap_jaccard |>
  dplyr::group_by(
    baseline_cluster
  ) |>
  dplyr::summarise(
    mean_Jaccard =
      mean(Jaccard),
    median_Jaccard =
      stats::median(Jaccard),
    min_Jaccard =
      min(Jaccard),
    q025_Jaccard =
      as.numeric(
        stats::quantile(
          Jaccard,
          0.025
        )
      ),
    q975_Jaccard =
      as.numeric(
        stats::quantile(
          Jaccard,
          0.975
        )
      ),
    max_Jaccard =
      max(Jaccard),
    .groups = "drop"
  )


cat("\n========================================\n")
cat("PARTICIPANT-LEVEL BOOTSTRAP STABILITY\n")
cat("========================================\n")

print(
  bootstrap_summary,
  n = Inf,
  width = Inf
)


cat("\n========================================\n")
cat("BOOTSTRAP CLUSTER-SPECIFIC JACCARD\n")
cat("========================================\n")

print(
  bootstrap_jaccard_summary,
  n = Inf,
  width = Inf
)


readr::write_csv(
  bootstrap_results,
  file.path(
    tables_dir,
    "cluster_bootstrap_ARI_1000.csv"
  )
)

readr::write_csv(
  bootstrap_summary,
  file.path(
    tables_dir,
    "cluster_bootstrap_summary.csv"
  )
)

readr::write_csv(
  bootstrap_jaccard,
  file.path(
    tables_dir,
    "cluster_bootstrap_Jaccard_1000.csv"
  )
)

readr::write_csv(
  bootstrap_jaccard_summary,
  file.path(
    tables_dir,
    "cluster_bootstrap_Jaccard_summary.csv"
  )
)


p_bootstrap_ari <- ggplot2::ggplot(
  bootstrap_results,
  ggplot2::aes(
    x = ARI
  )
) +
  ggplot2::geom_histogram(
    bins = 30
  ) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  ggplot2::labs(
    x = "Adjusted Rand Index vs baseline",
    y = "Bootstrap iterations",
    title = "Participant-level bootstrap stability"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "cluster_bootstrap_ARI_distribution.png"
  ),
  p_bootstrap_ari,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# G. CONSOLIDATED REPORTING OUTPUTS
# ============================================================
#
# The detailed CSV files above remain the primary numerical
# record. This compact long-format table collects the key
# values most likely to be used in the manuscript, supplement,
# response to reviewers, and Quarto reproducibility report.
# ============================================================

pc_exact_values <- pc_results$n_PC[
  pc_results$identical_to_baseline
]

consensus_between <- consensus_summary |>
  dplyr::filter(
    !baseline_same_cluster
  )

consensus_within <- consensus_summary |>
  dplyr::filter(
    baseline_same_cluster
  )


robustness_summary <- dplyr::bind_rows(

  baseline_internal_characterization |>
    dplyr::transmute(
      section = "Additional internal characterization (full SNV spectral space)",
      metric = metric,
      value_numeric = value,
      value_text = format(
        value,
        digits = 10
      )
    ),


  tibble::tibble(
    section = "Baseline",
    metric = c(
      "Retained PCs",
      "Cumulative variance explained (%)",
      "Baseline k",
      "Cluster sizes",
      "Average silhouette",
      "Participants with negative silhouette"
    ),
    value_numeric = c(
      N_PC_BASELINE,
      baseline_cumulative_variance,
      K_BASELINE,
      NA_real_,
      baseline_average_silhouette,
      sum(
        baseline_sil_df$silhouette_width < 0
      )
    ),
    value_text = c(
      as.character(N_PC_BASELINE),
      format(
        baseline_cumulative_variance,
        digits = 10
      ),
      as.character(K_BASELINE),
      paste(
        baseline_sizes,
        collapse = "/"
      ),
      format(
        baseline_average_silhouette,
        digits = 10
      ),
      as.character(
        sum(
          baseline_sil_df$silhouette_width < 0
        )
      )
    )
  ),

  tibble::tibble(
    section = "PCA dimensionality sensitivity",
    metric = c(
      "PC counts with exact baseline recovery",
      "Number of tested PC counts",
      "Number with exact baseline recovery"
    ),
    value_numeric = c(
      NA_real_,
      length(PC_VALUES),
      length(pc_exact_values)
    ),
    value_text = c(
      paste(
        pc_exact_values,
        collapse = ","
      ),
      as.character(
        length(PC_VALUES)
      ),
      as.character(
        length(pc_exact_values)
      )
    )
  ),

  tibble::tibble(
    section = "k = 3 versus k = 2",
    metric = c(
      "Silhouette k=2",
      "Silhouette k=3",
      "Calinski-Harabasz k=2",
      "Calinski-Harabasz k=3",
      "Davies-Bouldin k=2",
      "Davies-Bouldin k=3"
    ),
    value_numeric = c(
      k2_row$average_silhouette,
      k3_row$average_silhouette,
      k2_row$calinski_harabasz,
      k3_row$calinski_harabasz,
      k2_row$davies_bouldin,
      k3_row$davies_bouldin
    ),
    value_text = format(
      c(
        k2_row$average_silhouette,
        k3_row$average_silhouette,
        k2_row$calinski_harabasz,
        k3_row$calinski_harabasz,
        k2_row$davies_bouldin,
        k3_row$davies_bouldin
      ),
      digits = 10
    )
  ),

  tibble::tibble(
    section = "Leave-one-out",
    metric = c(
      "Mean ARI",
      "Median ARI",
      "Exact recovery proportion",
      "ARI >= 0.8 proportion"
    ),
    value_numeric = c(
      loo_summary$mean_ARI,
      loo_summary$median_ARI,
      loo_summary$proportion_ARI_1,
      loo_summary$proportion_ARI_ge_0_8
    ),
    value_text = format(
      c(
        loo_summary$mean_ARI,
        loo_summary$median_ARI,
        loo_summary$proportion_ARI_1,
        loo_summary$proportion_ARI_ge_0_8
      ),
      digits = 10
    )
  ),

  tibble::tibble(
    section = "80% subsampling",
    metric = c(
      "Mean ARI",
      "Median ARI",
      "Exact recovery proportion",
      "ARI >= 0.8 proportion"
    ),
    value_numeric = c(
      subsampling_summary$mean_ARI,
      subsampling_summary$median_ARI,
      subsampling_summary$proportion_ARI_1,
      subsampling_summary$proportion_ARI_ge_0_8
    ),
    value_text = format(
      c(
        subsampling_summary$mean_ARI,
        subsampling_summary$median_ARI,
        subsampling_summary$proportion_ARI_1,
        subsampling_summary$proportion_ARI_ge_0_8
      ),
      digits = 10
    )
  ),

  tibble::tibble(
    section = "Pairwise consensus",
    metric = c(
      "Within-cluster mean consensus",
      "Within-cluster median consensus",
      "Between-cluster mean consensus",
      "Between-cluster median consensus"
    ),
    value_numeric = c(
      consensus_within$mean_consensus,
      consensus_within$median_consensus,
      consensus_between$mean_consensus,
      consensus_between$median_consensus
    ),
    value_text = format(
      c(
        consensus_within$mean_consensus,
        consensus_within$median_consensus,
        consensus_between$mean_consensus,
        consensus_between$median_consensus
      ),
      digits = 10
    )
  ),

  tibble::tibble(
    section = "Bootstrap",
    metric = c(
      "Mean ARI",
      "Median ARI",
      "Exact recovery proportion",
      "ARI >= 0.8 proportion",
      "Mean unique participants per bootstrap sample"
    ),
    value_numeric = c(
      bootstrap_summary$mean_ARI,
      bootstrap_summary$median_ARI,
      bootstrap_summary$proportion_ARI_1,
      bootstrap_summary$proportion_ARI_ge_0_8,
      bootstrap_summary$mean_unique_participants
    ),
    value_text = format(
      c(
        bootstrap_summary$mean_ARI,
        bootstrap_summary$median_ARI,
        bootstrap_summary$proportion_ARI_1,
        bootstrap_summary$proportion_ARI_ge_0_8,
        bootstrap_summary$mean_unique_participants
      ),
      digits = 10
    )
  ),

  bootstrap_jaccard_summary |>
    dplyr::transmute(
      section = "Bootstrap cluster-specific Jaccard",
      metric = paste0(
        "Cluster ",
        baseline_cluster,
        " mean Jaccard"
      ),
      value_numeric = mean_Jaccard,
      value_text = format(
        mean_Jaccard,
        digits = 10
      )
    ),

  cluster_silhouette_summary |>
    dplyr::transmute(
      section = "Cluster-specific silhouette",
      metric = paste0(
        "Cluster ",
        cluster,
        " mean silhouette"
      ),
      value_numeric = mean_silhouette,
      value_text = format(
        mean_silhouette,
        digits = 10
      )
    )
)


readr::write_csv(
  robustness_summary,
  file.path(
    tables_dir,
    "cluster_robustness_summary.csv"
  )
)


# Save all compact result objects together for later Quarto use.

robustness_objects <- list(
  settings = list(
    seed = SEED,
    baseline_n_pc = N_PC_BASELINE,
    baseline_k = K_BASELINE,
    pc_values = PC_VALUES,
    k_values = K_VALUES,
    n_subsample = N_SUBSAMPLE,
    subsample_fraction = SUBSAMPLE_FRACTION,
    n_bootstrap = N_BOOTSTRAP
  ),
  pc_results = pc_results,
  k_results = k_results,
  k3_vs_k2 = k3_vs_k2,
  baseline_silhouette_by_participant =
    baseline_sil_df,
  baseline_silhouette_by_cluster =
    cluster_silhouette_summary,
  full_snv_silhouette_by_participant =
    full_snv_sil_df,
  baseline_internal_characterization =
    baseline_internal_characterization,
  permanova_full_snv =
    permanova_full_snv,
  betadisper_full_snv =
    betadisper_full_snv,
  betadisper_perm =
    betadisper_perm,
  loo_results = loo_results,
  loo_summary = loo_summary,
  loo_influence = loo_influence,
  subsampling_results =
    subsampling_results,
  subsampling_summary =
    subsampling_summary,
  consensus_matrix =
    consensus_matrix,
  consensus_summary =
    consensus_summary,
  bootstrap_results =
    bootstrap_results,
  bootstrap_summary =
    bootstrap_summary,
  bootstrap_jaccard =
    bootstrap_jaccard,
  bootstrap_jaccard_summary =
    bootstrap_jaccard_summary,
  robustness_summary =
    robustness_summary
)


saveRDS(
  robustness_objects,
  file.path(
    processed_dir,
    "cluster_robustness_objects.rds"
  )
)


# Individual reusable RDS outputs.

saveRDS(
  baseline_internal_characterization,
  file.path(
    processed_dir,
    "baseline_internal_characterization_full_snv.rds"
  )
)

saveRDS(
  pc_results,
  file.path(
    processed_dir,
    "cluster_PC_count_sensitivity.rds"
  )
)

saveRDS(
  k_results,
  file.path(
    processed_dir,
    "cluster_number_sensitivity_k2_k6.rds"
  )
)

saveRDS(
  loo_results,
  file.path(
    processed_dir,
    "cluster_LOO_results.rds"
  )
)

saveRDS(
  subsampling_results,
  file.path(
    processed_dir,
    "cluster_subsampling_ARI_1000.rds"
  )
)

saveRDS(
  bootstrap_results,
  file.path(
    processed_dir,
    "cluster_bootstrap_ARI_1000.rds"
  )
)


# ------------------------------------------------------------
# H. Reproducibility environment
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    tables_dir,
    "sessionInfo_cluster_stability.txt"
  )
)


# ------------------------------------------------------------
# I. Final console report
# ------------------------------------------------------------

cat("\n\n")
cat("========================================\n")
cat("CLUSTER ROBUSTNESS ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "Baseline: ",
  N_PC_BASELINE,
  " PCs, k = ",
  K_BASELINE,
  ", cluster sizes ",
  paste(
    baseline_sizes,
    collapse = " / "
  ),
  "\n",
  sep = ""
)

cat(
  "Baseline cumulative variance: ",
  round(
    baseline_cumulative_variance,
    3
  ),
  "%\n",
  sep = ""
)

cat(
  "Baseline average silhouette: ",
  round(
    baseline_average_silhouette,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "Negative baseline silhouette widths: ",
  sum(
    baseline_sil_df$silhouette_width < 0
  ),
  "\n",
  sep = ""
)

cat(
  "Full-SNV average silhouette: ",
  round(
    full_snv_average_silhouette,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "Full-SNV PERMANOVA R2 / F / p: ",
  round(
    permanova_R2,
    3
  ),
  " / ",
  round(
    permanova_F,
    2
  ),
  " / ",
  signif(
    permanova_p,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Full-SNV betadisper permutation p: ",
  signif(
    betadisper_perm_p,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "PC counts with exact baseline recovery: ",
  paste(
    pc_exact_values,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "LOO mean / median ARI: ",
  round(
    loo_summary$mean_ARI,
    3
  ),
  " / ",
  round(
    loo_summary$median_ARI,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "80% subsampling mean / median ARI: ",
  round(
    subsampling_summary$mean_ARI,
    3
  ),
  " / ",
  round(
    subsampling_summary$median_ARI,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "Bootstrap mean / median ARI: ",
  round(
    bootstrap_summary$mean_ARI,
    3
  ),
  " / ",
  round(
    bootstrap_summary$median_ARI,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "\nKey reporting table:\n",
  normalizePath(
    file.path(
      tables_dir,
      "cluster_robustness_summary.csv"
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
      "cluster_robustness_objects.rds"
    ),
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat("========================================\n")
