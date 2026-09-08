# ============================================================
# 08_tables_figures.R
#
# Applied Sciences revision
# FINAL publication tables / figures assembled from frozen 01-07 results.
#
# IMPORTANT
# ---------
# This script introduces no new inferential analyses.
# It assembles the final main and supplementary outputs from the frozen
# pediatric analyses and the final harmonized adult-reference branch.
#
# Main table:
#   Table 1  Pediatric cohort characteristics
#
# Main figures:
#   Figure 1  Curated analytical workflow asset (copied, not generated)
#   Figure 2  PCA / clustering / validation
#   Figure 3  Clinical composition
#   Figure 4  Supplementary references + clean-strip reproducibility
#   Figure 5  D1-D5 distributions
#   Figure 6  Raw absorbance spectra + IQR + D1-D5 regions
#   Figure 7  Signed standardized D1-D5 signatures
#
# Supplementary:
#   Table S1 adult domain-level clinical/spectral map
#   Table S2 cluster validation
#   Table S3 sex sensitivity
#   Table S4 strip technical reproducibility
#   Table S5/S5b substrate-reference invariance
#   Figure S1 consensus heatmap
#   Figure S2 sex-residualized PCA sensitivity
#   Figure S3 substrate / D1 sensitivity
#   Figure S4 adult PC2 / PC4 / PC5 loading profiles
#   Figure S5 strip pointwise technical variability
#
# Run from project root:
#   source("R/08_tables_figures.R")
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
  "ggplot2",
  "patchwork",
  "tibble",
  "dplyr",
  "tidyr",
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

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them and rerun R/08_tables_figures.R."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
})

set.seed(SEED)

processed_dir <- file.path(
  "data",
  "pediatric",
  "processed"
)

adult_processed_dir <- file.path(
  "data",
  "adult",
  "processed"
)

tables_dir <- "tables"
figures_dir <- "figures"
supplement_dir <- "supplement"

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

dir.create(
  supplement_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------

require_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required frozen result is missing: ",
      path,
      "\nRun the corresponding analysis script from 01-07 first."
    )
  }
  path
}

read_required_csv <- function(path) {
  readr::read_csv(
    require_file(path),
    show_col_types = FALSE
  )
}

save_figure <- function(
    plot,
    stem,
    width,
    height,
    dpi = 400
) {

  png_path <- file.path(
    figures_dir,
    paste0(stem, ".png")
  )

  pdf_path <- file.path(
    figures_dir,
    paste0(stem, ".pdf")
  )

  ggplot2::ggsave(
    png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    pdf_path,
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  c(
    png_path,
    pdf_path
  )
}

clean_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x
}

cluster_factor <- function(x) {
  factor(
    as.integer(as.character(x)),
    levels = 1:3,
    labels = c(
      "Cluster 1",
      "Cluster 2",
      "Cluster 3"
    )
  )
}

cluster_cols <- c(
  "Cluster 1" = "#C44E52",
  "Cluster 2" = "#4C72B0",
  "Cluster 3" = "#55A868"
)

# Supplementary-only palettes. These affect display only; no analytical
# objects, statistics, or numerical outputs are changed.
supp_sensitivity_cols <- c(
  "Baseline" = "#4C78A8",
  "Strip-aligned component removed" = "#F58518",
  "D1-dominant region excluded" = "#54A24B"
)

supp_pc_cols <- c(
  "PC2" = "#4C72B0",
  "PC4" = "#DD8452",
  "PC5" = "#55A868"
)

theme_pub <- function(base_size = 11) {
  ggplot2::theme_minimal(
    base_size = base_size
  ) +
    ggplot2::theme(
      panel.grid.minor =
        ggplot2::element_blank(),
      plot.title =
        ggplot2::element_text(
          face = "bold"
        ),
      strip.text =
        ggplot2::element_text(
          face = "bold"
        ),
      legend.title =
        ggplot2::element_text(
          face = "bold"
        )
    )
}


# ------------------------------------------------------------
# 2. Load frozen 01-07 objects
# ------------------------------------------------------------

X_snv <- as.matrix(
  readRDS(
    require_file(
      file.path(
        processed_dir,
        "pediatric_X_snv.rds"
      )
    )
  )
)

wn <- as.numeric(
  readRDS(
    require_file(
      file.path(
        processed_dir,
        "pediatric_wavenumbers.rds"
      )
    )
  )
)

cluster_baseline_raw <- readRDS(
  require_file(
    file.path(
      processed_dir,
      "pediatric_clusters_baseline.rds"
    )
  )
)

spectra_long <- readRDS(
  require_file(
    file.path(
      processed_dir,
      "pediatric_source_spectra_long.rds"
    )
  )
)

cluster_robustness <- readRDS(
  require_file(
    file.path(
      processed_dir,
      "cluster_robustness_objects.rds"
    )
  )
)

sex_obj <- readRDS(
  require_file(
    file.path(
      processed_dir,
      "sex_sensitivity_objects.rds"
    )
  )
)

strip_obj <- readRDS(
  require_file(
    file.path(
      processed_dir,
      "strip_reproducibility_objects.rds"
    )
  )
)

domain_obj <- readRDS(
  require_file(
    file.path(
      processed_dir,
      "spectral_domain_objects.rds"
    )
  )
)

# Figure 7 uses the frozen cluster-level standardized D1-D5 signatures
# produced by R/06_spectral_domains.R. Keep this as a reporting alias only;
# no values are recalculated here.
domain_signatures <- domain_obj$cluster_domain_signatures

if (is.null(domain_signatures)) {
  stop(
    "spectral_domain_objects.rds does not contain cluster_domain_signatures. ",
    "Run R/06_spectral_domains.R first."
  )
}

required_domain_signature_cols <- c(
  "cluster", "domain", "standardized_effect"
)

if (!all(required_domain_signature_cols %in% names(domain_signatures))) {
  stop("Frozen cluster-domain signature object has unexpected columns.")
}

# R/06 stores cluster as levels 1/2/3; convert once to the publication labels
# expected by the Figure 7 plotting code.
domain_signatures <- domain_signatures |>
  dplyr::mutate(
    cluster = cluster_factor(cluster)
  )

adult_obj <- readRDS(
  require_file(
    file.path(
      adult_processed_dir,
      "adult_reference_analysis_objects.rds"
    )
  )
)

if (!identical(dim(X_snv), c(25L, 793L))) {
  stop(
    "Expected pediatric X_snv = 25 x 793; obtained ",
    paste(dim(X_snv), collapse = " x "),
    "."
  )
}

if (length(wn) != ncol(X_snv)) {
  stop("Pediatric wavenumber vector does not match X_snv.")
}

sample_ids <- rownames(X_snv)

if (is.null(sample_ids)) {
  stop("X_snv must have participant IDs as row names.")
}

cluster_baseline <- cluster_baseline_raw

if (is.factor(cluster_baseline)) {
  cluster_baseline <- as.integer(
    as.character(cluster_baseline)
  )
} else {
  cluster_baseline <- as.integer(
    cluster_baseline
  )
}

if (!is.null(names(cluster_baseline_raw))) {
  cluster_baseline <- as.integer(
    cluster_baseline_raw[
      sample_ids
    ]
  )
}

names(cluster_baseline) <- sample_ids

if (!identical(
  as.integer(
    table(
      factor(
        cluster_baseline,
        levels = 1:3
      )
    )
  ),
  c(15L, 4L, 6L)
)) {
  stop(
    "Baseline cluster sizes are not 15 / 4 / 6."
  )
}


# ------------------------------------------------------------
# 3. Frozen reporting tables from 03-07
# ------------------------------------------------------------

pc_sensitivity <- read_required_csv(
  file.path(
    tables_dir,
    "cluster_PC_count_sensitivity.csv"
  )
)

k_sensitivity <- read_required_csv(
  file.path(
    tables_dir,
    "cluster_number_sensitivity_k2_k6.csv"
  )
)

baseline_silhouette <- read_required_csv(
  file.path(
    tables_dir,
    "baseline_silhouette_by_participant.csv"
  )
)

robustness_summary <- read_required_csv(
  file.path(
    tables_dir,
    "cluster_robustness_summary.csv"
  )
)

consensus_csv <- read_required_csv(
  file.path(
    tables_dir,
    "cluster_consensus_matrix_80pct.csv"
  )
)

sex_summary <- read_required_csv(
  file.path(
    tables_dir,
    "sex_sensitivity_summary.csv"
  )
)

adult_domain_map <- read_required_csv(
  file.path(
    tables_dir,
    "adult_domain_clinical_mapping.csv"
  )
)

adult_loading_aligned <- read_required_csv(
  file.path(
    tables_dir,
    "adult_PCA_variable_coordinates_aligned_PC2_PC4_PC5.csv"
  )
)

adult_axis_clinical_summary <- read_required_csv(
  file.path(
    tables_dir,
    "adult_PCA_selected_axis_clinical_summary.csv"
  )
)

adult_axis_category_means <- read_required_csv(
  file.path(
    tables_dir,
    "adult_PCA_selected_axis_category_means.csv"
  )
)

adult_axis_matches <- read_required_csv(
  file.path(
    tables_dir,
    "adult_PCA_selected_axis_matches.csv"
  )
)

adult_axis_score_correlations <- read_required_csv(
  file.path(
    tables_dir,
    "adult_PCA_historical_harmonized_score_correlations.csv"
  )
)

# Final adult-branch integrity checks. These do not perform new analyses;
# they only prevent stale or mismatched 07 outputs from being assembled.
if (!identical(as.integer(adult_obj$settings$n_reference_spectra), 18L)) {
  stop("Adult reference object does not report exactly 18 spectra.")
}

if (!identical(as.integer(adult_obj$settings$retained_variables), 793L)) {
  stop("Adult reference object does not report exactly 793 retained variables.")
}

required_adult_loading_cols <- c("wavenumber", "PC2", "PC4", "PC5")
if (!all(required_adult_loading_cols %in% names(adult_loading_aligned)) ||
    nrow(adult_loading_aligned) != 793L) {
  stop("Aligned adult loading export is not the expected 793-row PC2/PC4/PC5 table.")
}

expected_selected_axis_matches <- tibble::tribble(
  ~historical_PC, ~harmonized_PC,
  "PC1", "PC2",
  "PC2", "PC1",
  "PC4", "PC4",
  "PC5", "PC5"
)

observed_selected_axis_matches <- adult_axis_matches |>
  dplyr::select(historical_PC, harmonized_PC) |>
  dplyr::arrange(historical_PC)

if (!identical(
  as.data.frame(observed_selected_axis_matches),
  as.data.frame(expected_selected_axis_matches |> dplyr::arrange(historical_PC))
)) {
  stop("Adult historical-to-harmonized axis matching is not the expected PC1<->PC2 swap with PC4/PC5 preservation.")
}

expected_adult_map <- tibble::tribble(
  ~domain, ~PC, ~direction_class,
  "D1", "PC2", "predominantly negative",
  "D2", "PC2", "predominantly negative",
  "D3", "PC4", "predominantly positive",
  "D4", "PC4", "mixed",
  "D5", "PC5", "predominantly negative"
)

adult_map_check <- adult_domain_map |>
  dplyr::select(domain, PC, direction_class) |>
  dplyr::arrange(domain)

if (!identical(
  as.data.frame(adult_map_check),
  as.data.frame(expected_adult_map |> dplyr::arrange(domain))
)) {
  stop("Adult D1-D5 domain map does not match the final PC2/PC4/PC5 specification.")
}


# ------------------------------------------------------------
# 4. Reconstruct baseline PCA for display only
# ------------------------------------------------------------
#
# This uses the frozen X_snv matrix and the exact baseline specification.
# It is not a new analysis; it only regenerates plotting coordinates.

pca_baseline <- FactoMineR::PCA(
  X_snv,
  scale.unit = FALSE,
  ncp = 5,
  graph = FALSE
)

pca_scores <- tibble::as_tibble(
  as.data.frame(
    pca_baseline$ind$coord
  ),
  rownames = "participant_id"
) |>
  dplyr::mutate(
    cluster =
      cluster_factor(
        cluster_baseline[
          participant_id
        ]
      )
  )

pc1_var <- pca_baseline$eig[1, 2]
pc2_var <- pca_baseline$eig[2, 2]
pc3_var <- pca_baseline$eig[3, 2]

cum3 <- sum(
  pca_baseline$eig[
    1:3,
    2
  ]
)

if (abs(cum3 - 89.282) > 0.05) {
  stop(
    "Unexpected baseline PCA cumulative variance: ",
    round(cum3, 3),
    "%."
  )
}


# ============================================================
# MAIN TABLES
# ============================================================


# ------------------------------------------------------------
# 5. Table 1: individual pediatric cohort characteristics
# ------------------------------------------------------------

metadata <- sex_obj$metadata

required_metadata_cols <- c(
  "participant_id",
  "age",
  "sex",
  "pubertal_development",
  "endocrine_status",
  "uctd_status",
  "periodontal_status",
  "cluster"
)

if (!all(
  required_metadata_cols %in%
    names(metadata)
)) {
  stop(
    "sex_sensitivity_objects.rds does not contain the expected metadata columns."
  )
}

Table1 <- metadata |>
  dplyr::transmute(
    ID = participant_id,
    Age = age,
    Sex = as.character(sex),
    `Pubertal Development` =
      clean_label(
        pubertal_development
      ),
    `Endocrine Status` =
      clean_label(
        endocrine_status
      ),
    `UCTD Status` =
      clean_label(
        uctd_status
      ),
    `Periodontal Status` =
      clean_label(
        periodontal_status
      )
  ) |>
  dplyr::arrange(
    suppressWarnings(
      as.numeric(ID)
    ),
    ID
  )

readr::write_csv(
  Table1,
  file.path(
    tables_dir,
    "Table1_pediatric_cohort_characteristics.csv"
  )
)


# ============================================================
# SUPPLEMENTARY TABLES / REVIEWER-FACING EXPORTS
# ============================================================


# ------------------------------------------------------------
# 6. Table S1: adult domain-level clinical/spectral map
# ------------------------------------------------------------

readr::write_csv(
  adult_domain_map,
  file.path(
    supplement_dir,
    "TableS1_adult_domain_clinical_mapping.csv"
  )
)

# Aggregate adult axis-level clinical summaries are exported separately
# for reviewer inspection; they contain no participant identifiers.
readr::write_csv(
  adult_axis_clinical_summary,
  file.path(
    supplement_dir,
    "adult_PCA_selected_axis_clinical_summary.csv"
  )
)

readr::write_csv(
  adult_axis_category_means,
  file.path(
    supplement_dir,
    "adult_PCA_selected_axis_category_means.csv"
  )
)

readr::write_csv(
  adult_axis_matches,
  file.path(
    supplement_dir,
    "adult_PCA_selected_axis_matches.csv"
  )
)

readr::write_csv(
  adult_axis_score_correlations,
  file.path(
    supplement_dir,
    "adult_PCA_historical_harmonized_score_correlations.csv"
  )
)


# ------------------------------------------------------------
# 7. Table S2: cluster validation summary
# ------------------------------------------------------------

readr::write_csv(
  robustness_summary,
  file.path(
    supplement_dir,
    "TableS2_cluster_validation_summary.csv"
  )
)


# ------------------------------------------------------------
# 8. Table S3: sex sensitivity summary
# ------------------------------------------------------------

readr::write_csv(
  sex_summary,
  file.path(
    supplement_dir,
    "TableS3_sex_sensitivity_summary.csv"
  )
)


# ------------------------------------------------------------
# 9. Table S4: strip technical reproducibility
# ------------------------------------------------------------

strip_technical_summary <- strip_obj$technical_summary

readr::write_csv(
  strip_technical_summary,
  file.path(
    supplement_dir,
    "TableS4_strip_technical_reproducibility.csv"
  )
)


# ------------------------------------------------------------
# 10. Table S5 / S5b: strip-reference invariance
# ------------------------------------------------------------

strip_reference_invariance <- strip_obj$reference_residualization_summary

readr::write_csv(
  strip_reference_invariance,
  file.path(
    supplement_dir,
    "TableS5_strip_reference_invariance.csv"
  )
)

readr::write_csv(
  strip_obj$reference_residualization_pairwise_ARI,
  file.path(
    supplement_dir,
    "TableS5b_strip_reference_residualization_pairwise_ARI.csv"
  )
)


# ------------------------------------------------------------
# 11. Adult loading exports requested by reviewer
# ------------------------------------------------------------

adult_raw_loading_source <- require_file(
  file.path(
    tables_dir,
    "adult_PCA_variable_coordinates_raw_PC1_PC5.csv"
  )
)

adult_aligned_loading_source <- require_file(
  file.path(
    tables_dir,
    "adult_PCA_variable_coordinates_aligned_PC2_PC4_PC5.csv"
  )
)

file.copy(
  adult_raw_loading_source,
  file.path(
    supplement_dir,
    "adult_PCA_raw_variable_coordinates_PC1_PC5.csv"
  ),
  overwrite = TRUE
)

file.copy(
  adult_aligned_loading_source,
  file.path(
    supplement_dir,
    "adult_PCA_aligned_variable_coordinates_PC2_PC4_PC5.csv"
  ),
  overwrite = TRUE
)


# ============================================================
# FIGURE 1. CURATED ANALYTICAL WORKFLOW ASSET
# ============================================================

# Figure 1 is curated externally and is intentionally NOT regenerated here.
# The final asset must already contain the corrected analytical windows,
# PCA-based hierarchical clustering terminology, internal validation,
# and exploratory adult-reference domain mapping.

fig1_png <- file.path(
  figures_dir,
  "Figure1_analytical_workflow.png"
)

fig1_source <- file.path(
  figures_dir,
  "source_assets",
  "Figure1_analytical_workflow.png"
)

if (!file.exists(fig1_source)) {
  stop(
    "Curated Figure 1 is missing. Place the final asset at ",
    "figures/source_assets/Figure1_analytical_workflow.png before running 08."
  )
}

if (!isTRUE(file.copy(
  fig1_source,
  fig1_png,
  overwrite = TRUE
))) {
  stop("Failed to copy curated Figure 1 into figures/.")
}

fig1_paths <- fig1_png


# ============================================================
# FIGURE 2. PCA / CLUSTERING / VALIDATION
# ============================================================


# ------------------------------------------------------------
# 14. Panel A: baseline PCA
# ------------------------------------------------------------

p2_A <- ggplot2::ggplot(
  pca_scores,
  ggplot2::aes(
    x = Dim.1,
    y = Dim.2,
    color = cluster
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    color = "grey85"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    color = "grey85"
  ) +
  ggplot2::geom_point(
    size = 3.2,
    alpha = 0.90
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(
      group = cluster
    ),
    type = "norm",
    linewidth = 0.7,
    level = 0.70,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = cluster_cols
  ) +
  ggplot2::labs(
    title = "A. Baseline PCA / three-cluster solution",
    x = paste0(
      "PC1 (",
      round(pc1_var, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pc2_var, 1),
      "%)"
    ),
    color = "Cluster"
  ) +
  theme_pub()


# ------------------------------------------------------------
# 15. Panel B: k = 2-6 internal criteria
# ------------------------------------------------------------

k_long <- k_sensitivity |>
  dplyr::select(
    k,
    average_silhouette,
    calinski_harabasz,
    davies_bouldin
  ) |>
  tidyr::pivot_longer(
    cols = -k,
    names_to = "criterion",
    values_to = "value"
  ) |>
  dplyr::mutate(
    criterion =
      dplyr::recode(
        criterion,
        average_silhouette =
          "Average silhouette",
        calinski_harabasz =
          "Calinski-Harabasz",
        davies_bouldin =
          "Davies-Bouldin"
      )
  )

p2_B <- ggplot2::ggplot(
  k_long,
  ggplot2::aes(
    x = k,
    y = value
  )
) +
  ggplot2::geom_vline(
    xintercept = 3,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::geom_line(
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    size = 2.2
  ) +
  ggplot2::facet_wrap(
    ~ criterion,
    scales = "free_y",
    ncol = 1
  ) +
  ggplot2::scale_x_continuous(
    breaks = 2:6
  ) +
  ggplot2::labs(
    title = "B. Cluster-number sensitivity",
    x = "Number of clusters (k)",
    y = NULL
  ) +
  theme_pub(
    base_size = 9.5
  )


# ------------------------------------------------------------
# 16. Panel C: retained-PC sensitivity
# ------------------------------------------------------------

p2_C <- ggplot2::ggplot(
  pc_sensitivity,
  ggplot2::aes(
    x = n_PC,
    y = ARI_vs_baseline_3PC
  )
) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = 2,
    linewidth = 0.5
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 2.6
  ) +
  ggplot2::scale_x_continuous(
    breaks = 2:7
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      0,
      1.03
    )
  ) +
  ggplot2::labs(
    title = "C. Sensitivity to retained PCA dimensionality",
    x = "Number of retained PCs",
    y = "ARI vs baseline 3-PC partition"
  ) +
  theme_pub()


# ------------------------------------------------------------
# 17. Panel D: baseline silhouette
# ------------------------------------------------------------

baseline_silhouette <- baseline_silhouette |>
  dplyr::mutate(
    cluster =
      cluster_factor(cluster)
  )

p2_D <- ggplot2::ggplot(
  baseline_silhouette,
  ggplot2::aes(
    x = cluster,
    y = silhouette_width,
    color = cluster
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    linetype = 2
  ) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    size = 2.1
  ) +
  ggplot2::scale_color_manual(
    values = cluster_cols,
    guide = "none"
  ) +
  ggplot2::labs(
    title = "D. Baseline cluster-wise silhouette",
    x = NULL,
    y = "Silhouette width"
  ) +
  theme_pub()

Figure2 <- (
  p2_A |
    p2_B
) /
  (
    p2_C |
      p2_D
  ) +
  patchwork::plot_layout(
    widths = c(
      1.15,
      0.85
    )
  )

fig2_paths <- save_figure(
  Figure2,
  "Figure2_PCA_cluster_validation",
  width = 13,
  height = 10
)


# ============================================================
# FIGURE 3. CLINICAL COMPOSITION
# ============================================================


# ------------------------------------------------------------
# 18. Prepare clinical plotting data
# ------------------------------------------------------------

clinical_df <- metadata |>
  dplyr::mutate(
    cluster =
      cluster_factor(cluster),
    sex =
      factor(
        as.character(sex),
        levels = c(
          "F",
          "M"
        ),
        labels = c(
          "Female",
          "Male"
        )
      )
  )

p3_A <- ggplot2::ggplot(
  clinical_df,
  ggplot2::aes(
    x = cluster,
    fill = sex
  )
) +
  ggplot2::geom_bar(
    width = 0.70
  ) +
  ggplot2::labs(
    title = "A. Sex composition",
    x = NULL,
    y = "Participants",
    fill = "Sex"
  ) +
  theme_pub()

p3_B <- ggplot2::ggplot(
  clinical_df,
  ggplot2::aes(
    x = cluster,
    y = age,
    color = cluster
  )
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    size = 2.1
  ) +
  ggplot2::scale_color_manual(
    values = cluster_cols,
    guide = "none"
  ) +
  ggplot2::labs(
    title = "B. Age distribution",
    x = NULL,
    y = "Age, years"
  ) +
  theme_pub()

clinical_status_long <- clinical_df |>
  dplyr::select(
    participant_id,
    cluster,
    endocrine_status,
    uctd_status,
    periodontal_status
  ) |>
  tidyr::pivot_longer(
    cols = c(
      endocrine_status,
      uctd_status,
      periodontal_status
    ),
    names_to = "variable",
    values_to = "status"
  ) |>
  dplyr::mutate(
    variable =
      dplyr::recode(
        variable,
        endocrine_status =
          "Endocrine status",
        uctd_status =
          "UCTD status",
        periodontal_status =
          "Periodontal status"
      ),
    status =
      clean_label(status)
  )

p3_C <- ggplot2::ggplot(
  clinical_status_long,
  ggplot2::aes(
    x = cluster,
    fill = status
  )
) +
  ggplot2::geom_bar(
    width = 0.70
  ) +
  ggplot2::facet_wrap(
    ~ variable,
    nrow = 1
  ) +
  ggplot2::labs(
    title = "C-E. Documented clinical-status composition",
    x = NULL,
    y = "Participants",
    fill = "Status"
  ) +
  theme_pub(
    base_size = 9.5
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 20,
        hjust = 1
      )
  )

Figure3 <- (
  p3_A |
    p3_B
) /
  p3_C +
  patchwork::plot_layout(
    heights = c(
      0.85,
      1.15
    )
  )

fig3_paths <- save_figure(
  Figure3,
  "Figure3_clinical_composition",
  width = 13,
  height = 9
)


# ============================================================
# FIGURE 4. SUPPLEMENTARY REFERENCES / STRIP REPRODUCIBILITY
# ============================================================


# ------------------------------------------------------------
# 19. Reconstruct supplementary projection on frozen PCA axes
# ------------------------------------------------------------

supp_ids <- c(
  "etalon",
  "G",
  "N-1",
  "N-2",
  "N"
)

extract_source_spectrum <- function(id) {

  z <- spectra_long |>
    dplyr::filter(
      sample_id == id,
      (
        wavenumber >= 2800 &
          wavenumber <= 3400
      ) |
        (
          wavenumber >= 870 &
            wavenumber <= 1800
        )
    ) |>
    dplyr::arrange(
      dplyr::desc(
        wavenumber
      )
    )

  if (nrow(z) != length(wn)) {
    stop(
      "Unexpected retained spectral grid for supplementary sample ",
      id,
      ": ",
      nrow(z),
      " points."
    )
  }

  if (!isTRUE(
    all.equal(
      as.numeric(
        z$wavenumber
      ),
      wn,
      tolerance = 1e-6
    )
  )) {
    stop(
      "Supplementary wavenumber grid mismatch for ",
      id,
      "."
    )
  }

  as.numeric(
    z$absorbance
  )
}

X_sup_raw <- do.call(
  rbind,
  lapply(
    supp_ids,
    extract_source_spectrum
  )
)

rownames(X_sup_raw) <- supp_ids
colnames(X_sup_raw) <- colnames(X_snv)

X_sup_snv <- snv_matrix(
  X_sup_raw
)

X_strips_snv <- strip_obj$
  X_strips_snv

if (is.null(
  rownames(
    X_strips_snv
  )
)) {
  rownames(
    X_strips_snv
  ) <- paste0(
    "strip_",
    sprintf(
      "%02d",
      seq_len(
        nrow(
          X_strips_snv
        )
      )
    )
  )
}

if (!identical(
  dim(X_strips_snv),
  c(4L, 793L)
)) {
  stop(
    "Expected four independent strip spectra on the 793-variable grid."
  )
}

X_projection <- rbind(
  X_snv,
  X_sup_snv,
  X_strips_snv
)

ind_sup <- (
  nrow(X_snv) + 1L
):nrow(X_projection)

pca_projection <- FactoMineR::PCA(
  X_projection,
  scale.unit = FALSE,
  ncp = 5,
  ind.sup = ind_sup,
  graph = FALSE
)

active_projection <- tibble::as_tibble(
  as.data.frame(
    pca_projection$ind$coord
  ),
  rownames = "sample_id"
) |>
  dplyr::mutate(
    cluster =
      cluster_factor(
        cluster_baseline[
          sample_id
        ]
      ),
    role =
      "Active pediatric GCF"
  )

supp_projection <- tibble::as_tibble(
  as.data.frame(
    pca_projection$ind.sup$coord
  ),
  rownames = "sample_id"
) |>
  dplyr::mutate(
    reference_type =
      dplyr::case_when(
        sample_id %in%
          c(
            "N-1",
            "N-2",
            "N"
          ) ~
          "No documented pathology reference",
        sample_id == "G" ~
          "Gingivitis reference",
        sample_id == "etalon" ~
          "Archived clean strip",
        grepl(
          "^strip_",
          sample_id
        ) ~
          "Independent clean strip",
        TRUE ~
          "Supplementary reference"
      ),
    label =
      dplyr::case_when(
        sample_id == "etalon" ~
          "Archived strip",
        sample_id == "G" ~
          "G",
        TRUE ~
          sample_id
      )
  )


# ------------------------------------------------------------
# 20. Figure 4A: PCA projection
# ------------------------------------------------------------

p4_A <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    color = "grey85"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    color = "grey85"
  ) +
  ggplot2::geom_point(
    data = active_projection,
    ggplot2::aes(
      x = Dim.1,
      y = Dim.2,
      color = cluster
    ),
    size = 2.8,
    alpha = 0.80
  ) +
  ggplot2::geom_point(
    data = supp_projection,
    ggplot2::aes(
      x = Dim.1,
      y = Dim.2,
      shape = reference_type
    ),
    size = 4,
    color = "black",
    fill = "white",
    stroke = 1.0
  ) +
  ggplot2::geom_text(
    data = supp_projection,
    ggplot2::aes(
      x = Dim.1,
      y = Dim.2,
      label = label
    ),
    nudge_y = 0.6,
    size = 2.7,
    check_overlap = TRUE
  ) +
  ggplot2::scale_color_manual(
    values = cluster_cols
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "No documented pathology reference" = 21,
      "Gingivitis reference" = 24,
      "Archived clean strip" = 22,
      "Independent clean strip" = 23,
      "Supplementary reference" = 25
    )
  ) +
  ggplot2::labs(
    title = "A. Active pediatric spectra and supplementary references",
    x = paste0(
      "PC1 (",
      round(
        pca_projection$eig[1, 2],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        pca_projection$eig[2, 2],
        1
      ),
      "%)"
    ),
    color = "Cluster",
    shape = "Reference"
  ) +
  theme_pub(
    base_size = 10
  )


# ------------------------------------------------------------
# 21. Figure 4B: every clean reference vs cluster centroids
# ------------------------------------------------------------

strip_distance_df <- strip_obj$
  strip_reference_cluster_distances |>
  dplyr::mutate(
    cluster =
      factor(
        cluster,
        levels = 1:3,
        labels = c(
          "Cluster 1",
          "Cluster 2",
          "Cluster 3"
        )
      ),
    reference =
      clean_label(reference)
  )

p4_B <- ggplot2::ggplot(
  strip_distance_df,
  ggplot2::aes(
    x = cluster,
    y = distance,
    group = reference
  )
) +
  ggplot2::geom_line(
    linewidth = 0.45,
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    size = 2.0,
    alpha = 0.80
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(
      group = 1
    ),
    fun = stats::median,
    geom = "point",
    shape = 23,
    size = 4,
    fill = "white",
    color = "black"
  ) +
  ggplot2::labs(
    title = "B. Distance of clean-strip references to cluster centroids",
    subtitle = "All individual clean references identify Cluster 1 as nearest",
    x = NULL,
    y = "Euclidean distance in first five pediatric PCs"
  ) +
  theme_pub(
    base_size = 9.5
  )


# ------------------------------------------------------------
# 22. Figure 4C: four unused strips + archived reference
# ------------------------------------------------------------

X_strips_raw <- strip_obj$
  X_strips

if (!identical(
  dim(X_strips_raw),
  c(4L, 793L)
)) {
  stop(
    "Expected strip_obj$X_strips = 4 x 793."
  )
}

if (is.null(
  rownames(
    X_strips_raw
  )
)) {
  rownames(
    X_strips_raw
  ) <- rownames(
    X_strips_snv
  )
}

archived_strip_raw <- extract_source_spectrum(
  "etalon"
)

strip_overlay_matrix <- rbind(
  X_strips_raw,
  archived_etalon =
    archived_strip_raw
)

strip_overlay <- tibble::tibble(
  reference =
    rep(
      rownames(
        strip_overlay_matrix
      ),
      each =
        ncol(
          strip_overlay_matrix
        )
    ),
  wavenumber =
    rep(
      wn,
      times =
        nrow(
          strip_overlay_matrix
        )
    ),
  absorbance =
    as.vector(
      t(
        strip_overlay_matrix
      )
    )
) |>
  dplyr::mutate(
    window =
      dplyr::if_else(
        wavenumber > 2000,
        "3400-2800 cm^-1",
        "1800-870 cm^-1"
      )
  )

metric_value <- function(metric_name) {
  strip_technical_summary$value[
    strip_technical_summary$metric ==
      metric_name
  ][1]
}

strip_caption <- paste0(
  "Four unused strips: minimum pairwise Pearson r = ",
  format(
    metric_value(
      "Minimum pairwise Pearson"
    ),
    digits = 7
  ),
  "; minimum Spearman rho = ",
  format(
    metric_value(
      "Minimum pairwise Spearman"
    ),
    digits = 7
  ),
  "; integrated absorbance CV = ",
  round(
    metric_value(
      "Integrated absorbance CV (%)"
    ),
    2
  ),
  "%; median pointwise CV = ",
  round(
    metric_value(
      "Median pointwise CV (%)"
    ),
    2
  ),
  "%."
)

p4_C <- ggplot2::ggplot(
  strip_overlay,
  ggplot2::aes(
    x = wavenumber,
    y = absorbance,
    group = reference,
    linetype = reference
  )
) +
  ggplot2::geom_line(
    linewidth = 0.55,
    alpha = 0.85
  ) +
  ggplot2::facet_wrap(
    ~ window,
    scales = "free_x",
    ncol = 1
  ) +
  ggplot2::scale_x_reverse() +
  ggplot2::labs(
    title = "C. Reproducibility of independent unused strips",
    x = expression(
      "Wavenumber (cm"^{-1}*")"
    ),
    y = "Absorbance",
    linetype = "Reference",
    caption = strip_caption
  ) +
  theme_pub(
    base_size = 9
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.caption =
      ggplot2::element_text(
        hjust = 0,
        size = 8.5
      )
  )

Figure4 <- p4_A /
  (
    p4_B |
      p4_C
  ) +
  patchwork::plot_layout(
    heights = c(
      1.0,
      1.15
    ),
    widths = c(
      0.85,
      1.15
    )
  )

fig4_paths <- save_figure(
  Figure4,
  "Figure4_reference_projection_strip_reproducibility",
  width = 14,
  height = 11
)


# ============================================================
# FIGURE 5. D1-D5 DISTRIBUTIONS
# ============================================================


# ------------------------------------------------------------
# 23. Domain-level participant distributions
# ------------------------------------------------------------

domain_long <- domain_obj$
  domain_long |>
  dplyr::mutate(
    cluster =
      cluster_factor(cluster),
    domain =
      factor(
        domain,
        levels = paste0(
          "D",
          1:5
        )
      )
  )

domain_kw <- domain_obj$
  domain_kw |>
  dplyr::mutate(
    domain =
      factor(
        domain,
        levels = paste0(
          "D",
          1:5
        )
      ),
    p_label =
      paste0(
        "Kruskal-Wallis p = ",
        format.pval(
          p,
          digits = 3,
          eps = 1e-4
        )
      )
  )

p5 <- ggplot2::ggplot(
  domain_long,
  ggplot2::aes(
    x = cluster,
    y = value,
    color = cluster
  )
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    size = 1.8,
    alpha = 0.80
  ) +
  ggplot2::facet_wrap(
    ~ domain,
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::scale_color_manual(
    values = cluster_cols,
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Discriminative spectral-domain measures across the three FTIR-derived phenotypes",
    subtitle = paste(
      paste0(
        as.character(
          domain_kw$domain
        ),
        ": ",
        domain_kw$p_label
      ),
      collapse = "   |   "
    ),
    x = NULL,
    y = "Mean SNV-normalized spectral intensity within domain"
  ) +
  theme_pub(
    base_size = 9.5
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 8.5
      )
  )

fig5_paths <- save_figure(
  p5,
  "Figure5_D1_D5_domain_distributions",
  width = 14,
  height = 5.8
)


# ============================================================
# FIGURE 6. RAW ABSORBANCE MEDIAN + IQR
# ============================================================


# ------------------------------------------------------------
# 24. Raw cluster spectra and domain shading
# ------------------------------------------------------------

raw_summary <- domain_obj$
  raw_spectral_summary |>
  dplyr::mutate(
    cluster =
      cluster_factor(cluster),
    window =
      dplyr::if_else(
        wavenumber > 2000,
        "3400-2800 cm^-1",
        "1800-870 cm^-1"
      )
  )

domain_plot_defs <- tibble::tribble(
  ~domain, ~lower, ~upper,
  "D1",  876, 1127,
  "D2", 1150, 1313,
  "D3", 1382, 1413,
  "D4", 1477, 1706,
  "D5", 2942, 3092
) |>
  dplyr::mutate(
    window =
      dplyr::if_else(
        lower > 2000,
        "3400-2800 cm^-1",
        "1800-870 cm^-1"
      )
  )

p6 <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = domain_plot_defs,
    ggplot2::aes(
      xmin = lower,
      xmax = upper,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = "grey60",
    alpha = 0.08
  ) +
  ggplot2::geom_ribbon(
    data = raw_summary,
    ggplot2::aes(
      x = wavenumber,
      ymin = q25_absorbance,
      ymax = q75_absorbance,
      fill = cluster,
      group = cluster
    ),
    alpha = 0.13,
    linewidth = 0
  ) +
  ggplot2::geom_line(
    data = raw_summary,
    ggplot2::aes(
      x = wavenumber,
      y = median_absorbance,
      color = cluster,
      group = cluster
    ),
    linewidth = 0.75
  ) +
  ggplot2::facet_wrap(
    ~ window,
    scales = "free_x",
    ncol = 1
  ) +
  ggplot2::scale_x_reverse() +
  ggplot2::scale_color_manual(
    values = cluster_cols
  ) +
  ggplot2::scale_fill_manual(
    values = cluster_cols
  ) +
  ggplot2::labs(
    title = NULL,
    subtitle = "Lines: cluster medians; ribbons: interquartile ranges; shaded intervals: D1-D5",
    x = expression(
      "Wavenumber (cm"^{-1}*")"
    ),
    y = "Absorbance",
    color = "Cluster",
    fill = "Cluster"
  ) +
  theme_pub(
    base_size = 10
  ) +
  ggplot2::theme(
    legend.position = "bottom"
  )

fig6_paths <- save_figure(
  p6,
  "Figure6_raw_median_IQR_spectra_domains",
  width = 12,
  height = 8
)


# ============================================================
# FIGURE 7. SIGNED STANDARDIZED DOMAIN SIGNATURES
# ============================================================


# ------------------------------------------------------------
# 25. Heatmap
# ------------------------------------------------------------

heat_df <- domain_signatures |>
  dplyr::mutate(
    domain =
      factor(
        domain,
        levels = paste0(
          "D",
          1:5
        )
      ),
    cluster =
      factor(
        cluster,
        levels = c(
          "Cluster 1",
          "Cluster 2",
          "Cluster 3"
        )
      )
  )

p7 <- ggplot2::ggplot(
  heat_df,
  ggplot2::aes(
    x = domain,
    y = cluster,
    fill = standardized_effect
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 1.0
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%+.2f",
        standardized_effect
      )
    ),
    size = 4.2
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  ggplot2::labs(
    title = "Signed standardized multidomain spectral signatures",
    x = "Spectral domain",
    y = NULL,
    fill = "Standardized\ndeviation"
  ) +
  theme_pub(
    base_size = 11
  ) +
  ggplot2::theme(
    panel.grid =
      ggplot2::element_blank()
  )

fig7_paths <- save_figure(
  p7,
  "Figure7_signed_domain_heatmap",
  width = 8,
  height = 4.6
)


# ============================================================
# SUPPLEMENTARY FIGURES
# ============================================================


# ------------------------------------------------------------
# 28. Figure S1: consensus heatmap
# ------------------------------------------------------------

consensus_matrix <- as.matrix(
  consensus_csv[
    ,
    -1,
    drop = FALSE
  ]
)

rownames(
  consensus_matrix
) <- consensus_csv[[1]]

storage.mode(
  consensus_matrix
) <- "numeric"

consensus_ids <- rownames(
  consensus_matrix
)

cluster_for_consensus <-
  cluster_baseline[
    consensus_ids
  ]

consensus_order <- order(
  cluster_for_consensus,
  consensus_ids
)

ordered_ids <-
  consensus_ids[
    consensus_order
  ]

consensus_ordered <-
  consensus_matrix[
    ordered_ids,
    ordered_ids,
    drop = FALSE
  ]

consensus_long <- as.data.frame(
  as.table(
    consensus_ordered
  )
)

names(
  consensus_long
) <- c(
  "sample_1",
  "sample_2",
  "consensus"
)

consensus_long$sample_1 <-
  factor(
    consensus_long$sample_1,
    levels = rev(
      ordered_ids
    )
  )

consensus_long$sample_2 <-
  factor(
    consensus_long$sample_2,
    levels = ordered_ids
  )

pS1 <- ggplot2::ggplot(
  consensus_long,
  ggplot2::aes(
    x = sample_2,
    y = sample_1,
    fill = consensus
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08306B",
    limits = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    title = "Repeated-subsampling consensus matrix",
    x = NULL,
    y = NULL,
    fill = "Consensus"
  ) +
  theme_pub(
    base_size = 8
  ) +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6
      ),
    axis.text.y =
      ggplot2::element_text(
        size = 6
      ),
    panel.grid =
      ggplot2::element_blank()
  )

figS1_paths <- save_figure(
  pS1,
  "FigureS1_consensus_heatmap",
  width = 8,
  height = 7.5
)


# ------------------------------------------------------------
# 29. Figure S2: sex-residualized clustering sensitivity
# ------------------------------------------------------------

sex_resid <- sex_obj$
  sex_residualized_assignments |>
  dplyr::mutate(
    baseline_cluster =
      cluster_factor(
        baseline_cluster
      ),
    sex =
      factor(
        sex,
        levels = c(
          "F",
          "M"
        ),
        labels = c(
          "Female",
          "Male"
        )
      )
  )

baseline_for_sex <- pca_scores |>
  dplyr::select(
    participant_id,
    Dim.1,
    Dim.2,
    cluster
  ) |>
  dplyr::left_join(
    clinical_df |>
      dplyr::select(
        participant_id,
        sex
      ),
    by = "participant_id"
  )

pS2_A <- ggplot2::ggplot(
  baseline_for_sex,
  ggplot2::aes(
    x = Dim.1,
    y = Dim.2,
    color = cluster,
    shape = sex
  )
) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::scale_color_manual(
    values = cluster_cols
  ) +
  ggplot2::labs(
    title = "A. Baseline PCA",
    x = "PC1",
    y = "PC2",
    color = "Baseline cluster",
    shape = "Sex"
  ) +
  theme_pub(
    base_size = 9.5
  )

pS2_B <- ggplot2::ggplot(
  sex_resid,
  ggplot2::aes(
    x = residualized_PC1,
    y = residualized_PC2,
    color = baseline_cluster,
    shape = sex
  )
) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::scale_color_manual(
    values = cluster_cols
  ) +
  ggplot2::labs(
    title = paste0(
      "B. PCA after removing sex-associated spectral component; ARI = ",
      sprintf(
        "%.3f",
        sex_obj$
          sex_residualized_summary$
          ARI_vs_baseline
      )
    ),
    x = "Residual PC1",
    y = "Residual PC2",
    color = "Baseline cluster",
    shape = "Sex"
  ) +
  theme_pub(
    base_size = 9.5
  )

FigureS2 <- pS2_A |
  pS2_B

figS2_paths <- save_figure(
  FigureS2,
  "FigureS2_sex_residualized_PCA_sensitivity",
  width = 12,
  height = 5.5
)


# ------------------------------------------------------------
# 30. Figure S3: substrate / D1 sensitivity
# ------------------------------------------------------------

substrate_sensitivity <- tibble::tibble(
  analysis = factor(
    c(
      "Baseline",
      "Strip-aligned component removed",
      "D1-dominant region excluded"
    ),
    levels = c(
      "Baseline",
      "Strip-aligned component removed",
      "D1-dominant region excluded"
    )
  ),
  ARI_vs_baseline = c(
    1,
    strip_obj$ARI_strip_resid,
    strip_obj$ARI_no_D1
  )
)

pS3 <- ggplot2::ggplot(
  substrate_sensitivity,
  ggplot2::aes(
    x = analysis,
    y = ARI_vs_baseline,
    fill = analysis
  )
) +
  ggplot2::geom_col(
    width = 0.62,
    color = "grey25",
    linewidth = 0.25
  ) +
  ggplot2::scale_fill_manual(
    values = supp_sensitivity_cols,
    guide = "none"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.3f",
        ARI_vs_baseline
      )
    ),
    vjust = -0.35,
    size = 4
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      0,
      1.08
    )
  ) +
  ggplot2::labs(
    title = "Sensitivity of the pediatric partition to substrate-associated signal",
    x = NULL,
    y = "ARI vs baseline partition"
  ) +
  theme_pub() +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 20,
        hjust = 1
      )
  )

figS3_paths <- save_figure(
  pS3,
  "FigureS3_substrate_D1_sensitivity",
  width = 8.5,
  height = 5.5
)


# ------------------------------------------------------------
# 29. Figure S4: adult aligned PC2 / PC4 / PC5 loading profiles
# ------------------------------------------------------------

adult_loading_long <- adult_loading_aligned |>
  tidyr::pivot_longer(
    cols = c(
      PC2,
      PC4,
      PC5
    ),
    names_to = "PC",
    values_to = "loading"
  ) |>
  dplyr::mutate(
    window =
      dplyr::if_else(
        wavenumber > 2000,
        "3400-2800 cm^-1",
        "1800-870 cm^-1"
      )
  )

adult_domains_for_plot <- domain_plot_defs |>
  dplyr::mutate(
    PC = dplyr::case_when(
      domain %in% c("D1", "D2") ~ "PC2",
      domain %in% c("D3", "D4") ~ "PC4",
      domain == "D5" ~ "PC5",
      TRUE ~ NA_character_
    ),
    window =
      dplyr::if_else(
        lower > 2000,
        "3400-2800 cm^-1",
        "1800-870 cm^-1"
      )
  )

pS4 <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = adult_domains_for_plot,
    ggplot2::aes(
      xmin = lower,
      xmax = upper,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = "grey60",
    alpha = 0.07
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  ggplot2::geom_line(
    data = adult_loading_long,
    ggplot2::aes(
      x = wavenumber,
      y = loading,
      color = PC
    ),
    linewidth = 0.75
  ) +
  ggplot2::scale_color_manual(
    values = supp_pc_cols,
    guide = "none"
  ) +
  ggplot2::facet_grid(
    PC ~ window,
    scales = "free_x"
  ) +
  ggplot2::scale_x_reverse() +
  ggplot2::labs(
    title = "Adult-reference PCA loading profiles used for domain-level cross-cohort interpretation",
    subtitle = "Harmonized windows; signs fixed to the stated clinical-pole convention; shaded intervals indicate pediatric D1-D5",
    x = expression(
      "Wavenumber (cm"^{-1}*")"
    ),
    y = "FactoMineR variable coordinate"
  ) +
  theme_pub(
    base_size = 9
  )

figS4_paths <- save_figure(
  pS4,
  "FigureS4_adult_PC2_PC4_PC5_loading_profiles",
  width = 12,
  height = 8.5
)


# ------------------------------------------------------------
# 32. Figure S5: pointwise technical variability
# ------------------------------------------------------------

pointwise_strip <- strip_obj$
  pointwise

pS5 <- ggplot2::ggplot(
  pointwise_strip,
  ggplot2::aes(
    x = wavenumber,
    y = cv_percent
  )
) +
  ggplot2::geom_line(
    linewidth = 0.75,
    color = "#4C72B0"
  ) +
  ggplot2::facet_wrap(
    ~ dplyr::if_else(
      wavenumber > 2000,
      "3400-2800 cm^-1",
      "1800-870 cm^-1"
    ),
    scales = "free_x",
    ncol = 1
  ) +
  ggplot2::scale_x_reverse() +
  ggplot2::labs(
    title = "Pointwise technical variability across four independent unused strips",
    x = expression(
      "Wavenumber (cm"^{-1}*")"
    ),
    y = "Coefficient of variation (%)"
  ) +
  theme_pub(
    base_size = 9.5
  )

figS5_paths <- save_figure(
  pS5,
  "FigureS5_strip_pointwise_CV",
  width = 9,
  height = 7
)


# ============================================================
# OUTPUT MANIFEST
# ============================================================

main_table_paths <- c(
  file.path(tables_dir, "Table1_pediatric_cohort_characteristics.csv")
)

supplementary_table_paths <- c(
  file.path(supplement_dir, "TableS1_adult_domain_clinical_mapping.csv"),
  file.path(supplement_dir, "TableS2_cluster_validation_summary.csv"),
  file.path(supplement_dir, "TableS3_sex_sensitivity_summary.csv"),
  file.path(supplement_dir, "TableS4_strip_technical_reproducibility.csv"),
  file.path(supplement_dir, "TableS5_strip_reference_invariance.csv"),
  file.path(supplement_dir, "TableS5b_strip_reference_residualization_pairwise_ARI.csv"),
  file.path(supplement_dir, "adult_PCA_selected_axis_clinical_summary.csv"),
  file.path(supplement_dir, "adult_PCA_selected_axis_category_means.csv"),
  file.path(supplement_dir, "adult_PCA_selected_axis_matches.csv"),
  file.path(supplement_dir, "adult_PCA_historical_harmonized_score_correlations.csv"),
  file.path(supplement_dir, "adult_PCA_raw_variable_coordinates_PC1_PC5.csv"),
  file.path(supplement_dir, "adult_PCA_aligned_variable_coordinates_PC2_PC4_PC5.csv")
)

main_figure_paths <- c(
  fig1_paths,
  fig2_paths,
  fig3_paths,
  fig4_paths,
  fig5_paths,
  fig6_paths,
  fig7_paths
)

supplementary_figure_paths <- c(
  figS1_paths,
  figS2_paths,
  figS3_paths,
  figS4_paths,
  figS5_paths
)

manifest <- dplyr::bind_rows(
  tibble::tibble(output_type = "main_table", path = main_table_paths),
  tibble::tibble(output_type = "supplementary_table", path = supplementary_table_paths),
  tibble::tibble(output_type = "main_figure", path = main_figure_paths),
  tibble::tibble(output_type = "supplementary_figure", path = supplementary_figure_paths)
)

missing_manifest_outputs <- manifest$path[!file.exists(manifest$path)]
if (length(missing_manifest_outputs) > 0L) {
  stop(
    "Final output manifest contains missing files:\n",
    paste(missing_manifest_outputs, collapse = "\n")
  )
}

readr::write_csv(
  manifest,
  file.path(tables_dir, "08_output_manifest.csv")
)

cat("\n")
cat("========================================\n")
cat("FINAL TABLES / FIGURES COMPLETE\n")
cat("========================================\n")
cat("Main table: Table 1 created.\n")
cat("Main figures: Figure 1-7 assembled/created.\n")
cat("Supplementary: Table S1-S5/S5b + adult aggregate/loading exports.\n")
cat("Supplementary figures: Figure S1-S5 created.\n")
cat("\nChecks:\n")
cat(
  "  Baseline cluster sizes: ",
  paste(
    as.integer(table(factor(cluster_baseline, levels = 1:3))),
    collapse = "/"
  ),
  "\n",
  sep = ""
)
cat(
  "  Baseline 3-PC cumulative variance: ",
  round(cum3, 3),
  "%\n",
  sep = ""
)
cat(
  "  Adult reference retained variables: ",
  adult_obj$settings$retained_variables,
  "\n",
  sep = ""
)
cat("  Adult axis matching: historical PC1->PC2; PC2->PC1; PC4->PC4; PC5->PC5.\n")
cat("  Adult domain map: D1/D2->PC2; D3/D4->PC4; D5->PC5.\n")
cat("\nManifest: tables/08_output_manifest.csv\n")
cat("========================================\n")
