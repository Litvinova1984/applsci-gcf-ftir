# ============================================================
# 02_pca_clustering.R
#
# Pediatric GCF FTIR
# Baseline PCA and hierarchical clustering
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
  "ggplot2",
  "tibble",
  "dplyr",
  "readr"
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


set.seed(SEED)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

x_file <- file.path(
  "data",
  "pediatric",
  "processed",
  "pediatric_X_snv.rds"
)

processed_dir <- file.path(
  "data",
  "pediatric",
  "processed"
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
    "pediatric_X_snv.rds not found. Run R/01_validate_public_inputs.R and verify the public processed inputs."
  )
}


# ------------------------------------------------------------
# 2. Load SNV-normalized spectra
# ------------------------------------------------------------

X_snv <- as.matrix(
  readRDS(x_file)
)

if (!identical(
  dim(X_snv),
  c(25L, 793L)
)) {
  stop(
    "Expected X_snv dimensions 25 x 793."
  )
}

if (any(!is.finite(X_snv))) {
  stop(
    "Non-finite values detected in X_snv."
  )
}


# ------------------------------------------------------------
# 3. Full PCA for scree/eigenvalue assessment
# ------------------------------------------------------------
#
# No feature-wise scaling is applied after SNV.
#
# A full PCA is calculated here only to evaluate the
# eigenvalue spectrum and cumulative explained variance.

max_pc <- min(
  nrow(X_snv) - 1L,
  ncol(X_snv)
)

pca_full <- FactoMineR::PCA(
  X_snv,
  scale.unit = FALSE,
  ncp = max_pc,
  graph = FALSE
)


eig_raw <- as.data.frame(
  pca_full$eig
)

eig_table <- tibble::tibble(
  PC = seq_len(nrow(eig_raw)),
  eigenvalue = eig_raw[[1]],
  variance_percent = eig_raw[[2]],
  cumulative_percent = eig_raw[[3]]
)

readr::write_csv(
  eig_table,
  file.path(
    tables_dir,
    "PCA_eigenvalues_full.csv"
  )
)


# Baseline number of principal components used throughout this script.
N_PC_BASELINE <- 3L


# ------------------------------------------------------------
# 4. Scree plot
# ------------------------------------------------------------

p_scree <- eig_table |>
  dplyr::filter(PC <= 12) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = PC,
      y = variance_percent
    )
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::geom_vline(
    xintercept = N_PC_BASELINE,
    linetype = "dashed"
  ) +
  ggplot2::scale_x_continuous(
    breaks = 1:12
  ) +
  ggplot2::labs(
    x = "Principal component",
    y = "Explained variance (%)",
    title = "Scree plot of SNV-normalized pediatric GCF FTIR spectra"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  )

ggplot2::ggsave(
  file.path(
    figures_dir,
    "PCA_scree_plot.png"
  ),
  p_scree,
  width = 8,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# ------------------------------------------------------------
# 5. Principal-component retention
# ------------------------------------------------------------
#
# Three principal components were retained for the baseline
# clustering analysis based on visual inspection of the scree
# plot together with cumulative explained variance.
#
# PC1-PC3 jointly explain approximately 89.3% of the total
# variance.
#
# Sensitivity of the clustering solution to alternative
# numbers of retained PCs is evaluated separately in
# 03_cluster_stability.R.


retention_summary <- eig_table |>
  dplyr::filter(PC <= 7)

readr::write_csv(
  retention_summary,
  file.path(
    tables_dir,
    "PCA_component_retention.csv"
  )
)


readr::write_csv(
  retention_summary,
  file.path(
    tables_dir,
    "PCA_component_retention.csv"
  )
)


# ------------------------------------------------------------
# ------------------------------------------------------------
# 6. Baseline PCA used for clustering
# ------------------------------------------------------------
#
# Three principal components are retained for clustering.
# No additional feature-wise scaling is applied.


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

# ------------------------------------------------------------
# 7. Baseline hierarchical clustering
# ------------------------------------------------------------
#
# The submitted baseline solution comprised three clusters.
# Here k = 3 is fixed to reproduce that baseline solution.
#
# Alternative cluster numbers (k = 2:6) are evaluated
# independently in the cluster-robustness analysis
# (03_cluster_stability.R).
#
# Clustering:
#   Euclidean distance
#   Ward.D2 linkage
#   first three PCA scores

K_BASELINE <- 3L

dist_baseline <- stats::dist(
  scores_baseline,
  method = "euclidean"
)

hc_baseline <- stats::hclust(
  dist_baseline,
  method = "ward.D2"
)

cluster_baseline <- stats::cutree(
  hc_baseline,
  k = K_BASELINE
)


# ------------------------------------------------------------
# 8. Baseline cluster membership
# ------------------------------------------------------------

cluster_table <- tibble::tibble(
  sample_id = names(cluster_baseline),
  cluster = as.integer(cluster_baseline)
)


readr::write_csv(
  cluster_table,
  file.path(
    tables_dir,
    "baseline_cluster_membership.csv"
  )
)


# ------------------------------------------------------------
# 9. Baseline cluster sizes
# ------------------------------------------------------------

cluster_sizes <- tibble::tibble(
  cluster = as.integer(
    names(
      table(cluster_baseline)
    )
  ),
  n = as.integer(
    table(cluster_baseline)
  )
)


readr::write_csv(
  cluster_sizes,
  file.path(
    tables_dir,
    "baseline_cluster_sizes.csv"
  )
)


# ------------------------------------------------------------
# 10. Save reusable analytical objects
# ------------------------------------------------------------

saveRDS(
  pca_baseline,
  file.path(
    processed_dir,
    "pediatric_PCA_baseline_3PC.rds"
  )
)

saveRDS(
  scores_baseline,
  file.path(
    processed_dir,
    "pediatric_PCA_scores_3PC.rds"
  )
)

# Baseline cluster assignments are consumed by scripts 03-06 and 08.
saveRDS(
  cluster_baseline,
  file.path(
    processed_dir,
    "pediatric_clusters_baseline.rds"
  )
)

# ------------------------------------------------------------
# 11. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    tables_dir,
    "sessionInfo_PCA_clustering.txt"
  )
)


# ------------------------------------------------------------
# 12. Console report
# ------------------------------------------------------------

cat("\n========================================\n")
cat("BASELINE PCA / CLUSTERING COMPLETE\n")
cat("========================================\n")


cat(
  "Input matrix: ",
  nrow(X_snv),
  " x ",
  ncol(X_snv),
  "\n",
  sep = ""
)


cat(
  "Additional feature scaling: NO\n"
)


cat(
  "Retained PCs: ",
  N_PC_BASELINE,
  "\n",
  sep = ""
)


cat(
  "Cumulative variance: ",
  round(
    eig_table$cumulative_percent[
      N_PC_BASELINE
    ],
    3
  ),
  "%\n",
  sep = ""
)


cat(
  "Distance: Euclidean\n"
)


cat(
  "Linkage: Ward.D2\n"
)


cat(
  "Baseline k: ",
  K_BASELINE,
  "\n",
  sep = ""
)


cat(
  "Cluster sizes: ",
  paste(
    as.integer(
      table(cluster_baseline)
    ),
    collapse = " / "
  ),
  "\n",
  sep = ""
)


cat("========================================\n")
