# ============================================================
# 07_tables_figures.R
#
# Pediatric GCF FTIR exploratory reframe
# Publication tables and figures from the pediatric-only workflow
#
# IMPORTANT
# ---------
# This script introduces no new model-selection or inferential branch.
# It assembles publication outputs from the frozen objects produced by
# scripts 01-06. The retained scientific scope is limited to:
#
#   - direct pediatric spectral variation and PCA profiles;
#   - an exploratory three-group partition;
#   - internal stability and specification sensitivity;
#   - sex-associated variation;
#   - clean-strip technical reproducibility and sensitivity;
#   - participant-level changes across sensitivity partitions.
#
# Run from the project root after scripts 01-06:
#   source("R/07_tables_figures.R")
# ============================================================


# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "00_utils.R"))

required_pkgs <- c(
  "ggplot2",
  "patchwork",
  "tibble",
  "dplyr",
  "tidyr",
  "readr"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them and rerun R/07_tables_figures.R."
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

processed_dir <- file.path("data", "pediatric", "processed")
tables_dir <- "tables"
figures_dir <- "figures"
supplement_dir <- "supplement"

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------

require_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required result is missing: ",
      path,
      "\nRun the corresponding script from 01-06 first."
    )
  }
  path
}

read_required_csv <- function(path) {
  readr::read_csv(require_file(path), show_col_types = FALSE)
}

save_figure <- function(plot, stem, width, height, dpi = 400) {
  png_path <- file.path(figures_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figures_dir, paste0(stem, ".pdf"))

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
    bg = "white"
  )

  c(png_path, pdf_path)
}

fmt_num <- function(x, digits = 3) {
  formatC(as.numeric(x), digits = digits, format = "f")
}

fmt_p <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(format.pval(x, digits = 2, eps = 0.001))
  formatC(x, digits = 3, format = "f")
}

format_n_pct <- function(n, total) {
  sprintf("%d (%.1f%%)", as.integer(n), 100 * as.numeric(n) / total)
}

group_factor <- function(x) {
  factor(
    as.integer(as.character(x)),
    levels = 1:3,
    labels = c("Group 1", "Group 2", "Group 3")
  )
}

partition_factor <- function(x) {
  factor(
    x,
    levels = c(
      "Baseline",
      "Sex-residualized",
      "Clean-strip-residualized"
    )
  )
}

spectral_window <- function(x) {
  factor(
    ifelse(as.numeric(x) >= 2800, "3400-2800", "1800-870"),
    levels = c("3400-2800", "1800-870")
  )
}

clean_reference_label <- function(x) {
  dplyr::recode(
    as.character(x),
    archived_clean_strip = "Archived strip",
    strip_01 = "Independent strip 1",
    strip_02 = "Independent strip 2",
    strip_03 = "Independent strip 3",
    strip_04 = "Independent strip 4",
    mean_four_strips = "Mean of four strips",
    .default = as.character(x)
  )
}

theme_pub <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.tag = ggplot2::element_text(face = "bold")
    )
}

group_cols <- c(
  "Group 1" = "#C44E52",
  "Group 2" = "#4C72B0",
  "Group 3" = "#55A868"
)

sex_cols <- c(
  "F" = "#B24C63",
  "M" = "#3F6FA0"
)

metric_lookup <- function(tbl, metric, value_col = "value") {
  z <- tbl[tbl$metric == metric, value_col, drop = TRUE]
  if (length(z) != 1L) {
    stop("Expected exactly one row for metric: ", metric)
  }
  as.numeric(z[[1]])
}

# Remove stale publication outputs from the previous analytical framing.
retired_paths <- c(
  file.path(tables_dir, "Table2_adult_informed_orientation_score.csv"),
  file.path(supplement_dir, "TableS1_adult_directional_map.csv"),
  file.path(supplement_dir, "adult_PCA_raw_loadings_PC1_PC4_PC5.csv"),
  file.path(supplement_dir, "adult_PCA_raw_loadings_PC1_PC5.csv")
)

retired_figure_stems <- c(
  "Figure2_PCA_cluster_validation",
  "Figure3_clinical_composition",
  "Figure4_reference_projection_strip_reproducibility",
  "Figure5_D1_D5_domain_distributions",
  "Figure6_raw_median_IQR_spectra_domains",
  "Figure7_signed_domain_heatmap",
  "Figure8_adult_informed_orientation_score",
  "Figure9_control_proximity_vs_adult_orientation",
  "FigureS2_sex_residualized_PCA_sensitivity",
  "FigureS3_substrate_D1_sensitivity",
  "FigureS4_adult_PC1_PC4_PC5_loading_profiles",
  "FigureS5_strip_pointwise_CV"
)

retired_paths <- c(
  retired_paths,
  unlist(
    lapply(
      retired_figure_stems,
      function(stem) file.path(figures_dir, paste0(stem, c(".png", ".pdf")))
    )
  )
)

unlink(retired_paths[file.exists(retired_paths)], force = TRUE)


# ------------------------------------------------------------
# 2. Load and validate frozen outputs from scripts 01-06
# ------------------------------------------------------------

X_raw <- as.matrix(
  readRDS(
    require_file(
      file.path(processed_dir, "pediatric_X_active_absorbance.rds")
    )
  )
)

X_snv <- as.matrix(
  readRDS(
    require_file(
      file.path(processed_dir, "pediatric_X_snv.rds")
    )
  )
)

wavenumbers <- as.numeric(
  readRDS(
    require_file(
      file.path(processed_dir, "pediatric_wavenumbers.rds")
    )
  )
)

pca_baseline <- readRDS(
  require_file(
    file.path(processed_dir, "pediatric_PCA_baseline_3PC.rds")
  )
)

scores_baseline <- as.matrix(
  readRDS(
    require_file(
      file.path(processed_dir, "pediatric_PCA_scores_3PC.rds")
    )
  )
)

pca_profiles <- readRDS(
  require_file(
    file.path(processed_dir, "pediatric_PCA_profiles_PC1_PC3.rds")
  )
)

cluster_baseline_raw <- readRDS(
  require_file(
    file.path(processed_dir, "pediatric_clusters_baseline.rds")
  )
)

robustness <- readRDS(
  require_file(
    file.path(processed_dir, "cluster_robustness_objects.rds")
  )
)

sex_obj <- readRDS(
  require_file(
    file.path(processed_dir, "sex_sensitivity_objects.rds")
  )
)

strip_obj <- readRDS(
  require_file(
    file.path(processed_dir, "strip_reproducibility_objects.rds")
  )
)

transition_obj <- readRDS(
  require_file(
    file.path(processed_dir, "assignment_transition_objects.rds")
  )
)

pca_sex_residualized <- readRDS(
  require_file(
    file.path(processed_dir, "pediatric_PCA_sex_residualized.rds")
  )
)

if (!identical(dim(X_raw), c(25L, 793L))) {
  stop("Expected raw pediatric matrix dimensions 25 x 793.")
}

if (!identical(dim(X_snv), c(25L, 793L))) {
  stop("Expected SNV pediatric matrix dimensions 25 x 793.")
}

if (!identical(rownames(X_raw), rownames(X_snv))) {
  stop("Raw and SNV matrices have different participant ordering.")
}

if (length(wavenumbers) != ncol(X_raw)) {
  stop("The wavenumber vector does not match the pediatric matrices.")
}

cluster_baseline <- if (is.factor(cluster_baseline_raw)) {
  as.integer(as.character(cluster_baseline_raw))
} else {
  as.integer(cluster_baseline_raw)
}

names(cluster_baseline) <- names(cluster_baseline_raw)

if (!identical(rownames(X_snv), names(cluster_baseline))) {
  stop("Baseline cluster assignments do not match matrix ordering.")
}

if (!identical(dim(scores_baseline), c(25L, 3L))) {
  stop("Expected baseline PCA score dimensions 25 x 3.")
}

if (
  length(pca_profiles$wavenumber) != 793L ||
    !identical(dim(pca_profiles$eigenvector_loadings), c(793L, 3L))
) {
  stop("Unexpected dimensions in the saved PC1-PC3 profile object.")
}

metadata <- sex_obj$metadata |>
  dplyr::mutate(
    participant_id = as.character(participant_id),
    baseline_cluster = as.integer(baseline_cluster)
  )

if (!setequal(metadata$participant_id, rownames(X_snv))) {
  stop("Metadata IDs do not match the active pediatric matrix.")
}

metadata <- metadata[
  match(rownames(X_snv), metadata$participant_id),
  ,
  drop = FALSE
]

n_total <- nrow(metadata)


# ============================================================
# MAIN TABLES
# ============================================================


# ------------------------------------------------------------
# 3. Table 1: pediatric cohort characteristics
# ------------------------------------------------------------

count_level <- function(variable, level) {
  sum(as.character(metadata[[variable]]) == level, na.rm = TRUE)
}

Table1 <- dplyr::bind_rows(
  tibble::tibble(
    Characteristic = c("Age, years", "Age range, years"),
    Category = c("Mean +/- SD", "Minimum-maximum"),
    Value = c(
      sprintf("%.2f +/- %.2f", mean(metadata$age), stats::sd(metadata$age)),
      paste0(min(metadata$age), "-", max(metadata$age))
    )
  ),
  tibble::tibble(
    Characteristic = "Sex",
    Category = c("Female", "Male"),
    Value = c(
      format_n_pct(count_level("sex", "F"), n_total),
      format_n_pct(count_level("sex", "M"), n_total)
    )
  ),
  tibble::tibble(
    Characteristic = "Pubertal development",
    Category = c("Normal", "Early", "Delayed", "Unknown / insufficient information"),
    Value = vapply(
      c("Normal", "Early", "Delayed", "Unknown"),
      function(z) format_n_pct(count_level("pubertal_development", z), n_total),
      character(1)
    )
  ),
  tibble::tibble(
    Characteristic = "Endocrine status",
    Category = c("Documented pathology", "No documented pathology", "Unknown"),
    Value = vapply(
      c("Documented pathology", "No documented pathology", "Unknown"),
      function(z) format_n_pct(count_level("endocrine_status", z), n_total),
      character(1)
    )
  ),
  tibble::tibble(
    Characteristic = "UCTD status",
    Category = c("Documented pathology", "No documented pathology"),
    Value = vapply(
      c("Documented pathology", "No documented pathology"),
      function(z) format_n_pct(count_level("uctd_status", z), n_total),
      character(1)
    )
  ),
  tibble::tibble(
    Characteristic = "Periodontal status",
    Category = c("Documented pathology", "No documented pathology"),
    Value = vapply(
      c("Documented pathology", "No documented pathology"),
      function(z) format_n_pct(count_level("periodontal_status", z), n_total),
      character(1)
    )
  )
)

readr::write_csv(
  Table1,
  file.path(tables_dir, "Table1_pediatric_cohort_characteristics.csv")
)


# ------------------------------------------------------------
# 4. Table 2: exploratory partition and internal stability
# ------------------------------------------------------------

baseline_variance <- as.data.frame(pca_baseline$eig)[3, 3]
baseline_sizes <- as.integer(table(factor(cluster_baseline, levels = 1:3)))

Table2 <- dplyr::bind_rows(
  tibble::tibble(
    Analysis = "Baseline representation",
    Metric = c(
      "Participants",
      "Spectral variables",
      "Retained principal components",
      "Cumulative variance explained (%)",
      "Exploratory group sizes",
      "Average silhouette"
    ),
    Value = c(
      as.character(nrow(X_snv)),
      as.character(ncol(X_snv)),
      as.character(robustness$settings$baseline_n_pc),
      fmt_num(baseline_variance, 2),
      paste(baseline_sizes, collapse = "/"),
      fmt_num(mean(robustness$baseline_silhouette_by_participant$silhouette_width), 3)
    )
  ),
  tibble::tibble(
    Analysis = "PCA dimensionality sensitivity",
    Metric = c(
      "Tested PC counts",
      "PC counts reproducing the baseline partition",
      "Number reproducing the baseline partition"
    ),
    Value = c(
      paste(robustness$pc_results$n_PC, collapse = ", "),
      paste(
        robustness$pc_results$n_PC[robustness$pc_results$identical_to_baseline],
        collapse = ", "
      ),
      as.character(sum(robustness$pc_results$identical_to_baseline))
    )
  ),
  tibble::tibble(
    Analysis = "Leave-one-out",
    Metric = c("Mean ARI", "Median ARI", "Exact recovery (%)"),
    Value = c(
      fmt_num(robustness$loo_summary$mean_ARI, 3),
      fmt_num(robustness$loo_summary$median_ARI, 3),
      fmt_num(100 * robustness$loo_summary$proportion_ARI_1, 1)
    )
  ),
  tibble::tibble(
    Analysis = "Repeated 80% subsampling",
    Metric = c(
      "Mean ARI",
      "Median ARI",
      "Within-group mean consensus",
      "Between-group mean consensus"
    ),
    Value = c(
      fmt_num(robustness$subsampling_summary$mean_ARI, 3),
      fmt_num(robustness$subsampling_summary$median_ARI, 3),
      fmt_num(
        robustness$consensus_summary$mean_consensus[
          robustness$consensus_summary$baseline_same_cluster
        ],
        3
      ),
      fmt_num(
        robustness$consensus_summary$mean_consensus[
          !robustness$consensus_summary$baseline_same_cluster
        ],
        3
      )
    )
  ),
  tibble::tibble(
    Analysis = "Participant bootstrap",
    Metric = c("Mean ARI", "Median ARI", "Exact recovery (%)"),
    Value = c(
      fmt_num(robustness$bootstrap_summary$mean_ARI, 3),
      fmt_num(robustness$bootstrap_summary$median_ARI, 3),
      fmt_num(100 * robustness$bootstrap_summary$proportion_ARI_1, 1)
    )
  ),
  robustness$bootstrap_jaccard_summary |>
    dplyr::transmute(
      Analysis = "Participant bootstrap",
      Metric = paste0("Group ", baseline_cluster, " mean Jaccard"),
      Value = fmt_num(mean_Jaccard, 3)
    )
)

readr::write_csv(
  Table2,
  file.path(tables_dir, "Table2_exploratory_partition_stability.csv")
)


# ------------------------------------------------------------
# 5. Table 3: sex-associated and clean-strip sensitivity
# ------------------------------------------------------------

pc1_row <- sex_obj$pc_sex_results |>
  dplyr::filter(PC == "PC1")

sex_res <- sex_obj$sex_residualized_summary
strip_technical <- strip_obj$technical_summary
strip_primary <- strip_obj$primary_sensitivity_summary
pairwise_partitions <- transition_obj$pairwise_partition_summary

Table3 <- dplyr::bind_rows(
  tibble::tibble(
    Source = "Sex-associated variation",
    Metric = c(
      "Sex x baseline partition Fisher p",
      "Sex x baseline partition Cramer's V",
      "Cramer's V 95% CI",
      "PC1 Wilcoxon raw p",
      "PC1 Wilcoxon BH-adjusted p",
      "PC1 rank-biserial correlation (F vs M)",
      "PC1 Hedges g (F minus M)"
    ),
    Value = c(
      fmt_p(sex_obj$sex_cluster_association$fisher_p_value),
      fmt_num(sex_obj$sex_cluster_association$cramers_v, 3),
      paste0(
        fmt_num(sex_obj$sex_cluster_association$cramers_v_ci_lower, 3),
        " to ",
        fmt_num(sex_obj$sex_cluster_association$cramers_v_ci_upper, 3)
      ),
      fmt_p(pc1_row$p_value_raw),
      fmt_p(pc1_row$p_value_BH),
      fmt_num(pc1_row$rank_biserial_F_vs_M, 3),
      fmt_num(pc1_row$hedges_g_F_minus_M, 3)
    )
  ),
  tibble::tibble(
    Source = "Sex-residualized sensitivity",
    Metric = c(
      "ARI vs baseline",
      "Aligned group sizes",
      "Participants changing aligned assignment",
      "Average silhouette",
      "Sex x residualized partition Fisher p"
    ),
    Value = c(
      fmt_num(sex_res$ARI_vs_baseline, 3),
      sex_res$residualized_cluster_sizes,
      as.character(sex_res$participants_changing_aligned_assignment),
      fmt_num(sex_res$residualized_average_silhouette, 3),
      fmt_p(sex_res$sex_cluster_fisher_p_after_residualization)
    )
  ),
  tibble::tibble(
    Source = "Clean-strip technical reproducibility",
    Metric = c(
      "Minimum pairwise Pearson r",
      "Minimum pairwise Spearman rho",
      "Integrated absorbance CV (%)",
      "Median pointwise CV (%)",
      "95th percentile pointwise CV (%)"
    ),
    Value = c(
      fmt_num(metric_lookup(strip_technical, "Minimum pairwise Pearson"), 6),
      fmt_num(metric_lookup(strip_technical, "Minimum pairwise Spearman"), 6),
      fmt_num(metric_lookup(strip_technical, "Integrated absorbance CV (%)"), 3),
      fmt_num(metric_lookup(strip_technical, "Median pointwise CV (%)"), 3),
      fmt_num(metric_lookup(strip_technical, "95th percentile pointwise CV (%)"), 3)
    )
  ),
  tibble::tibble(
    Source = "Clean-strip-direction sensitivity",
    Metric = c(
      "ARI vs baseline",
      "Participants changing aligned assignment",
      "Average silhouette"
    ),
    Value = c(
      fmt_num(metric_lookup(strip_primary, "ARI vs baseline"), 3),
      as.character(
        round(metric_lookup(strip_primary, "Participants changing aligned assignment"))
      ),
      fmt_num(metric_lookup(strip_primary, "Average silhouette"), 3)
    )
  ),
  tibble::tibble(
    Source = "Cross-sensitivity comparison",
    Metric = "ARI: sex-residualized vs clean-strip-residualized",
    Value = fmt_num(
      pairwise_partitions$ARI[
        pairwise_partitions$comparison ==
          "Sex-residualized vs Clean-strip-residualized"
      ],
      3
    )
  )
)

readr::write_csv(
  Table3,
  file.path(tables_dir, "Table3_sex_and_substrate_sensitivity.csv")
)


# ============================================================
# SUPPLEMENTARY TABLES
# ============================================================

readr::write_csv(
  robustness$robustness_summary,
  file.path(supplement_dir, "TableS1_cluster_stability_detailed.csv")
)

readr::write_csv(
  sex_obj$pc_sex_results,
  file.path(supplement_dir, "TableS2_PC1_PC3_sex_associations_BH.csv")
)

readr::write_csv(
  strip_obj$technical_summary,
  file.path(supplement_dir, "TableS3_clean_strip_technical_reproducibility.csv")
)

readr::write_csv(
  strip_obj$reference_residualization_summary,
  file.path(supplement_dir, "TableS4_clean_strip_direction_sensitivity.csv")
)

readr::write_csv(
  transition_obj$assignment_transitions,
  file.path(supplement_dir, "TableS5_participant_assignment_transitions.csv")
)

readr::write_csv(
  transition_obj$pairwise_partition_summary,
  file.path(supplement_dir, "TableS6_pairwise_partition_agreement.csv")
)

pca_loading_export <- tibble::tibble(
  wavenumber = pca_profiles$wavenumber,
  variable = colnames(X_snv),
  PC1 = pca_profiles$eigenvector_loadings[, 1],
  PC2 = pca_profiles$eigenvector_loadings[, 2],
  PC3 = pca_profiles$eigenvector_loadings[, 3]
) |>
  dplyr::arrange(dplyr::desc(wavenumber))

readr::write_csv(
  pca_loading_export,
  file.path(supplement_dir, "TableS7_PCA_eigenvector_loadings_PC1_PC3.csv")
)


# ============================================================
# FIGURE 1. ANALYTICAL WORKFLOW
# ============================================================

wrap_text <- function(x, width = 32) {
  paste(strwrap(x, width = width), collapse = "\n")
}

draw_box <- function(x, y, w, h, title, body, fill = "#F7F7F7") {
  grid::grid.roundrect(
    x = grid::unit(x, "npc"),
    y = grid::unit(y, "npc"),
    width = grid::unit(w, "npc"),
    height = grid::unit(h, "npc"),
    r = grid::unit(0.016, "npc"),
    gp = grid::gpar(fill = fill, col = "#4D4D4D", lwd = 1.3)
  )

  grid::grid.text(
    title,
    x = grid::unit(x, "npc"),
    y = grid::unit(y + h / 2 - 0.040, "npc"),
    gp = grid::gpar(fontsize = 10.5, fontface = "bold")
  )

  grid::grid.text(
    wrap_text(body, width = 31),
    x = grid::unit(x, "npc"),
    y = grid::unit(y - 0.018, "npc"),
    gp = grid::gpar(fontsize = 8.2, lineheight = 0.92)
  )
}

draw_arrow <- function(x1, y1, x2, y2) {
  grid::grid.segments(
    x0 = grid::unit(x1, "npc"),
    y0 = grid::unit(y1, "npc"),
    x1 = grid::unit(x2, "npc"),
    y1 = grid::unit(y2, "npc"),
    arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"),
    gp = grid::gpar(lwd = 1.4, col = "#4D4D4D")
  )
}

make_workflow_figure <- function() {
  grid::grid.newpage()

  grid::grid.text(
    "Exploratory analysis of pediatric GCF FTIR measurements",
    x = 0.5,
    y = 0.955,
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )

  grid::grid.text(
    "Pediatric-only workflow with explicit stability and technical-sensitivity assessment",
    x = 0.5,
    y = 0.915,
    gp = grid::gpar(fontsize = 10.5)
  )

  draw_box(
    0.14, 0.69, 0.22, 0.29,
    "1. Measurements",
    paste0(
      "Pediatric GCF collected on absorbent strips (n = 25). ",
      "FTIR windows 3400-2800 and 1800-870 cm^-1; 793 retained variables."
    ),
    "#EEF3F8"
  )

  draw_box(
    0.38, 0.69, 0.22, 0.29,
    "2. Preprocessing",
    paste0(
      "Participant-level spectra; three technical acquisitions averaged before analysis. ",
      "Row-wise SNV; no additional feature-wise scaling."
    ),
    "#F3EFF7"
  )

  draw_box(
    0.62, 0.69, 0.22, 0.29,
    "3. Exploratory structure",
    paste0(
      "PCA followed by Euclidean/Ward.D2 clustering in a fixed three-PC space. ",
      "The three-group partition is descriptive, not a biological classification."
    ),
    "#F3F6EE"
  )

  draw_box(
    0.86, 0.69, 0.22, 0.29,
    "4. Internal stability",
    paste0(
      "Alternative k and PC counts, leave-one-out, repeated 80% subsampling, ",
      "consensus analysis, participant bootstrap, and group-specific Jaccard stability."
    ),
    "#F8F3EA"
  )

  draw_arrow(0.25, 0.69, 0.27, 0.69)
  draw_arrow(0.49, 0.69, 0.51, 0.69)
  draw_arrow(0.73, 0.69, 0.75, 0.69)

  draw_box(
    0.31, 0.315, 0.42, 0.25,
    "Sex-associated sensitivity",
    paste0(
      "PC1-PC3 comparisons with BH adjustment and effect sizes. ",
      "Feature-wise sex residualization followed by a fresh PCA and the same fixed clustering procedure."
    ),
    "#F7F7F7"
  )

  draw_box(
    0.69, 0.315, 0.42, 0.25,
    "Clean-strip technical sensitivity",
    paste0(
      "Four independent unused strips plus the archived reference. ",
      "Technical reproducibility, supplementary projection, and removal of an independently defined strip-aligned direction."
    ),
    "#F7F7F7"
  )

  draw_arrow(0.52, 0.545, 0.34, 0.44)
  draw_arrow(0.76, 0.545, 0.68, 0.44)

  grid::grid.roundrect(
    x = 0.5,
    y = 0.085,
    width = 0.88,
    height = 0.105,
    r = grid::unit(0.015, "npc"),
    gp = grid::gpar(fill = "#FAFAFA", col = "#5A5A5A", lwd = 1.2)
  )

  grid::grid.text(
    paste0(
      "Outputs: continuous PC1-PC3 loading profiles; an exploratory partition with quantified stability;\n",
      "evidence that sex-associated and substrate-aligned sources of variation affect the observed organization."
    ),
    x = 0.5,
    y = 0.085,
    gp = grid::gpar(fontsize = 9.6, fontface = "bold")
  )
}

fig1_png <- file.path(figures_dir, "Figure1_analytical_workflow.png")
fig1_pdf <- file.path(figures_dir, "Figure1_analytical_workflow.pdf")

grDevices::png(fig1_png, width = 5600, height = 3600, res = 500)
make_workflow_figure()
grDevices::dev.off()

grDevices::pdf(fig1_pdf, width = 14, height = 9, useDingbats = FALSE)
make_workflow_figure()
grDevices::dev.off()

fig1_paths <- c(fig1_png, fig1_pdf)


# ============================================================
# FIGURE 2. SPECTRAL VARIATION, PCA SCORES, AND CONTINUOUS LOADINGS
# ============================================================

raw_long <- as.data.frame(X_raw, check.names = FALSE) |>
  tibble::rownames_to_column(var = "participant_id") |>
  tidyr::pivot_longer(
    cols = -participant_id,
    names_to = "variable",
    values_to = "absorbance"
  ) |>
  dplyr::mutate(
    wavenumber = extract_wavenumber(variable),
    window = spectral_window(wavenumber)
  )

raw_median <- raw_long |>
  dplyr::group_by(window, wavenumber) |>
  dplyr::summarise(
    median_absorbance = stats::median(absorbance),
    .groups = "drop"
  )

p2_raw <- ggplot2::ggplot(
  raw_long,
  ggplot2::aes(wavenumber, absorbance, group = participant_id)
) +
  ggplot2::geom_line(alpha = 0.18, linewidth = 0.25) +
  ggplot2::geom_line(
    data = raw_median,
    ggplot2::aes(wavenumber, median_absorbance, group = 1),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  ggplot2::scale_x_reverse() +
  ggplot2::facet_grid(
    . ~ window,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::labs(
    title = "Participant-level absorbance spectra",
    subtitle = "Thin lines: individual measurements; thick line: cohort median",
    x = expression("Wavenumber (cm"^{-1} * ")"),
    y = "Absorbance"
  ) +
  theme_pub(10)

score_df <- tibble::as_tibble(
  scores_baseline,
  rownames = "participant_id"
) |>
  dplyr::mutate(
    group = group_factor(cluster_baseline[participant_id])
  )

pc_var <- as.data.frame(pca_baseline$eig)

p2_scores <- ggplot2::ggplot(
  score_df,
  ggplot2::aes(PC1, PC2, color = group)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(size = 2.6, alpha = 0.9) +
  ggplot2::scale_color_manual(values = group_cols, drop = FALSE) +
  ggplot2::labs(
    title = "Exploratory three-group representation",
    x = paste0("PC1 (", fmt_num(pc_var[1, 2], 2), "%)"),
    y = paste0("PC2 (", fmt_num(pc_var[2, 2], 2), "%)"),
    color = "Exploratory group"
  ) +
  theme_pub(10)

loading_long <- pca_loading_export |>
  tidyr::pivot_longer(
    cols = c(PC1, PC2, PC3),
    names_to = "PC",
    values_to = "loading"
  ) |>
  dplyr::mutate(
    PC = factor(PC, levels = c("PC1", "PC2", "PC3")),
    window = spectral_window(wavenumber)
  )

p2_loadings <- ggplot2::ggplot(
  loading_long,
  ggplot2::aes(wavenumber, loading)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::scale_x_reverse() +
  ggplot2::facet_grid(
    PC ~ window,
    scales = "free",
    space = "free_x"
  ) +
  ggplot2::labs(
    title = "Continuous PC1-PC3 eigenvector loading profiles",
    subtitle = "Raw mathematical PCA signs are retained; no post hoc orientation is applied",
    x = expression("Wavenumber (cm"^{-1} * ")"),
    y = "Loading coefficient"
  ) +
  theme_pub(10)

Figure2 <- (
  p2_raw | p2_scores
) / p2_loadings +
  patchwork::plot_layout(heights = c(1, 1.45)) +
  patchwork::plot_annotation(tag_levels = "A")

fig2_paths <- save_figure(
  Figure2,
  "Figure2_spectral_structure_and_PCA_profiles",
  width = 14,
  height = 11
)


# ============================================================
# FIGURE 3. PARTITION SELECTION AND INTERNAL STABILITY
# ============================================================

k_plot_data <- robustness$k_results |>
  dplyr::select(
    k,
    average_silhouette,
    calinski_harabasz,
    davies_bouldin
  ) |>
  tidyr::pivot_longer(
    cols = -k,
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c(
        "average_silhouette",
        "calinski_harabasz",
        "davies_bouldin"
      ),
      labels = c(
        "Average silhouette",
        "Calinski-Harabasz",
        "Davies-Bouldin"
      )
    )
  )

p3_k <- ggplot2::ggplot(
  k_plot_data,
  ggplot2::aes(k, value)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_vline(xintercept = 3, linetype = "dashed") +
  ggplot2::scale_x_continuous(breaks = 2:6) +
  ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  ggplot2::labs(
    title = "Internal indices across candidate group counts",
    x = "Number of groups (k)",
    y = NULL
  ) +
  theme_pub(9)

p3_pc <- ggplot2::ggplot(
  robustness$pc_results,
  ggplot2::aes(n_PC, ARI_vs_baseline_3PC)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::geom_vline(xintercept = 3, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
  ggplot2::scale_x_continuous(breaks = 2:7) +
  ggplot2::scale_y_continuous(limits = c(0, 1.03)) +
  ggplot2::labs(
    title = "Sensitivity to the number of retained PCs",
    x = "Retained PCs",
    y = "ARI vs 3-PC baseline"
  ) +
  theme_pub(10)

silhouette_df <- robustness$baseline_silhouette_by_participant
if (!"sample_id" %in% names(silhouette_df)) {
  silhouette_df$sample_id <- rownames(X_snv)
}
silhouette_df <- silhouette_df |>
  dplyr::arrange(cluster, silhouette_width) |>
  dplyr::mutate(
    participant_order = factor(sample_id, levels = sample_id),
    group = group_factor(cluster)
  )

p3_sil <- ggplot2::ggplot(
  silhouette_df,
  ggplot2::aes(participant_order, silhouette_width, fill = group)
) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
  ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Participant-level silhouette widths",
    x = "Participant",
    y = "Silhouette width",
    fill = "Exploratory group"
  ) +
  theme_pub(9) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.5))

resampling_df <- dplyr::bind_rows(
  robustness$loo_results |>
    dplyr::transmute(method = "Leave-one-out", ARI = ARI),
  robustness$subsampling_results |>
    dplyr::transmute(method = "80% subsampling", ARI = ARI),
  robustness$bootstrap_results |>
    dplyr::transmute(method = "Participant bootstrap", ARI = ARI)
) |>
  dplyr::mutate(
    method = factor(
      method,
      levels = c(
        "Leave-one-out",
        "80% subsampling",
        "Participant bootstrap"
      )
    )
  )

p3_resampling <- ggplot2::ggplot(
  resampling_df,
  ggplot2::aes(method, ARI)
) +
  ggplot2::geom_violin(trim = TRUE, linewidth = 0.4) +
  ggplot2::geom_boxplot(width = 0.13, outlier.shape = NA) +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 2.2) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
  ggplot2::scale_y_continuous(limits = c(-0.05, 1.05)) +
  ggplot2::labs(
    title = "Resampling agreement with the baseline partition",
    x = NULL,
    y = "Adjusted Rand Index"
  ) +
  theme_pub(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))

Figure3 <- (
  (p3_k / p3_pc) | (p3_sil / p3_resampling)
) +
  patchwork::plot_layout(widths = c(1.25, 1)) +
  patchwork::plot_annotation(tag_levels = "A")

fig3_paths <- save_figure(
  Figure3,
  "Figure3_partition_selection_and_stability",
  width = 14,
  height = 11
)


# ============================================================
# FIGURE 4. SEX-ASSOCIATED VARIATION AND SENSITIVITY
# ============================================================

scores_metadata <- score_df |>
  dplyr::left_join(
    metadata |>
      dplyr::select(participant_id, sex, baseline_cluster),
    by = "participant_id"
  )

p4_baseline <- ggplot2::ggplot(
  scores_metadata,
  ggplot2::aes(PC1, PC2, color = sex, shape = group)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(size = 2.8, alpha = 0.9) +
  ggplot2::scale_color_manual(values = sex_cols) +
  ggplot2::labs(
    title = "Baseline PCA scores by sex",
    x = "PC1 score",
    y = "PC2 score",
    color = "Sex",
    shape = "Exploratory group"
  ) +
  theme_pub(10)

pc_label_map <- setNames(
  paste0(
    sex_obj$pc_sex_results$PC,
    "\nBH p = ",
    vapply(sex_obj$pc_sex_results$p_value_BH, fmt_p, character(1))
  ),
  sex_obj$pc_sex_results$PC
)

pc_score_long <- scores_metadata |>
  dplyr::select(participant_id, sex, PC1, PC2, PC3) |>
  tidyr::pivot_longer(
    cols = c(PC1, PC2, PC3),
    names_to = "PC",
    values_to = "score"
  ) |>
  dplyr::mutate(PC = factor(PC, levels = c("PC1", "PC2", "PC3")))

p4_pc <- ggplot2::ggplot(
  pc_score_long,
  ggplot2::aes(sex, score, color = sex)
) +
  ggplot2::geom_boxplot(width = 0.5, outlier.shape = NA, color = "grey35") +
  ggplot2::geom_jitter(width = 0.08, size = 1.8, alpha = 0.85) +
  ggplot2::scale_color_manual(values = sex_cols) +
  ggplot2::facet_wrap(
    ~ PC,
    scales = "free_y",
    labeller = ggplot2::as_labeller(pc_label_map)
  ) +
  ggplot2::labs(
    title = "PC1-PC3 score distributions",
    x = "Sex",
    y = "PCA score",
    color = "Sex"
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "none")

sex_resid_scores <- as.matrix(
  pca_sex_residualized$ind$coord[, 1:3, drop = FALSE]
)

sex_resid_assignments <- sex_obj$sex_residualized_assignments |>
  dplyr::mutate(participant_id = as.character(participant_id))

sex_resid_df <- tibble::as_tibble(
  sex_resid_scores,
  rownames = "participant_id"
) |>
  dplyr::left_join(
    sex_resid_assignments |>
      dplyr::select(participant_id, sex, sex_residualized_cluster),
    by = "participant_id"
  ) |>
  dplyr::mutate(group = group_factor(sex_residualized_cluster))

p4_resid <- ggplot2::ggplot(
  sex_resid_df,
  ggplot2::aes(Dim.1, Dim.2, color = group, shape = sex)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
  ggplot2::geom_point(size = 2.8, alpha = 0.9) +
  ggplot2::scale_color_manual(values = group_cols, drop = FALSE) +
  ggplot2::labs(
    title = "Partition after feature-wise sex residualization",
    subtitle = paste0("ARI vs baseline = ", fmt_num(sex_res$ARI_vs_baseline, 3)),
    x = "Residualized PC1 score",
    y = "Residualized PC2 score",
    color = "Aligned group",
    shape = "Sex"
  ) +
  theme_pub(10)

composition_baseline <- metadata |>
  dplyr::transmute(
    partition = "Baseline",
    group = group_factor(baseline_cluster),
    sex = as.character(sex)
  )

composition_resid <- sex_resid_assignments |>
  dplyr::transmute(
    partition = "Sex-residualized",
    group = group_factor(sex_residualized_cluster),
    sex = as.character(sex)
  )

composition_df <- dplyr::bind_rows(
  composition_baseline,
  composition_resid
) |>
  dplyr::count(partition, group, sex, name = "n") |>
  dplyr::group_by(partition, group) |>
  dplyr::mutate(proportion = n / sum(n)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    partition = factor(partition, levels = c("Baseline", "Sex-residualized"))
  )

p4_composition <- ggplot2::ggplot(
  composition_df,
  ggplot2::aes(group, proportion, fill = sex)
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = sex_cols) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  ggplot2::facet_wrap(~ partition) +
  ggplot2::labs(
    title = "Sex composition before and after residualization",
    x = NULL,
    y = "Within-group proportion",
    fill = "Sex"
  ) +
  theme_pub(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

Figure4 <- (
  p4_baseline | p4_pc
) / (
  p4_resid | p4_composition
) +
  patchwork::plot_annotation(tag_levels = "A")

fig4_paths <- save_figure(
  Figure4,
  "Figure4_sex_associated_variation",
  width = 14,
  height = 10.5
)


# ============================================================
# FIGURE 5. CLEAN-STRIP REPRODUCIBILITY AND SENSITIVITY
# ============================================================

strip_overlay <- dplyr::bind_rows(
  tibble::tibble(
    reference = "Archived strip",
    wavenumber = wavenumbers,
    absorbance = strip_obj$archived_strip_raw
  ),
  dplyr::bind_rows(
    lapply(
      seq_len(nrow(strip_obj$X_strips)),
      function(i) {
        tibble::tibble(
          reference = paste0("Independent strip ", i),
          wavenumber = wavenumbers,
          absorbance = as.numeric(strip_obj$X_strips[i, ])
        )
      }
    )
  )
) |>
  dplyr::mutate(window = spectral_window(wavenumber))

p5_overlay <- ggplot2::ggplot(
  strip_overlay,
  ggplot2::aes(wavenumber, absorbance, color = reference)
) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.85) +
  ggplot2::scale_x_reverse() +
  ggplot2::facet_grid(
    . ~ window,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::labs(
    title = "Archived and independently measured unused strips",
    x = expression("Wavenumber (cm"^{-1} * ")"),
    y = "Absorbance",
    color = NULL
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "bottom")

pointwise_cv <- strip_obj$pointwise_variability |>
  dplyr::mutate(window = spectral_window(wavenumber))

p5_cv <- ggplot2::ggplot(
  pointwise_cv,
  ggplot2::aes(wavenumber, cv_percent)
) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::scale_x_reverse() +
  ggplot2::facet_grid(
    . ~ window,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::labs(
    title = "Pointwise variability across four independent strips",
    x = expression("Wavenumber (cm"^{-1} * ")"),
    y = "CV (%)"
  ) +
  theme_pub(9)

active_projection <- as.matrix(
  strip_obj$pca_clean_projection$ind$coord[, 1:2, drop = FALSE]
)
reference_projection <- as.matrix(
  strip_obj$pca_clean_projection$ind.sup$coord[, 1:2, drop = FALSE]
)

projection_active_df <- tibble::as_tibble(
  active_projection,
  rownames = "participant_id"
) |>
  dplyr::mutate(
    group = group_factor(cluster_baseline[participant_id]),
    observation_type = "Pediatric measurement"
  )

projection_reference_df <- tibble::as_tibble(
  reference_projection,
  rownames = "reference"
) |>
  dplyr::mutate(
    reference_label = clean_reference_label(reference),
    observation_type = "Clean-strip reference"
  )

p5_projection <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = projection_active_df,
    ggplot2::aes(Dim.1, Dim.2, color = group),
    size = 2.2,
    alpha = 0.75
  ) +
  ggplot2::geom_point(
    data = projection_reference_df,
    ggplot2::aes(Dim.1, Dim.2),
    shape = 17,
    size = 3.0
  ) +
  ggplot2::geom_text(
    data = projection_reference_df,
    ggplot2::aes(Dim.1, Dim.2, label = reference_label),
    nudge_y = 0.35,
    size = 2.5,
    check_overlap = TRUE
  ) +
  ggplot2::scale_color_manual(values = group_cols, drop = FALSE) +
  ggplot2::labs(
    title = "Clean-strip references projected into the fixed pediatric PCA space",
    subtitle = "Reference spectra do not contribute to the PCA axes",
    x = "PC1 score",
    y = "PC2 score",
    color = "Exploratory group"
  ) +
  theme_pub(9)

strip_sensitivity_plot <- strip_obj$reference_residualization_summary |>
  dplyr::mutate(
    reference_label = clean_reference_label(reference),
    reference_label = factor(reference_label, levels = reference_label)
  )

p5_sensitivity <- ggplot2::ggplot(
  strip_sensitivity_plot,
  ggplot2::aes(reference_label, ARI_vs_baseline)
) +
  ggplot2::geom_segment(
    ggplot2::aes(xend = reference_label, y = 0, yend = ARI_vs_baseline),
    linewidth = 0.5
  ) +
  ggplot2::geom_point(size = 2.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0("changed n=", changed_participants)),
    vjust = -0.8,
    size = 2.7
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1.05)) +
  ggplot2::labs(
    title = "Sensitivity to the selected clean-strip direction",
    subtitle = "The same fixed 3-PC / Ward.D2 / k=3 procedure is repeated",
    x = NULL,
    y = "ARI vs baseline"
  ) +
  theme_pub(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

Figure5 <- (
  p5_overlay | p5_cv
) / (
  p5_projection | p5_sensitivity
) +
  patchwork::plot_annotation(tag_levels = "A")

fig5_paths <- save_figure(
  Figure5,
  "Figure5_clean_strip_reproducibility_and_sensitivity",
  width = 14,
  height = 10.5
)


# ============================================================
# FIGURE 6. PARTICIPANT-LEVEL PARTITION TRANSITIONS
# ============================================================

transition_wide <- transition_obj$assignment_transitions |>
  dplyr::mutate(
    participant_id = as.character(participant_id),
    id_numeric = suppressWarnings(as.numeric(participant_id))
  ) |>
  dplyr::arrange(baseline_cluster, sex, id_numeric, participant_id) |>
  dplyr::mutate(
    participant_label = paste0(participant_id, " (", sex, ")"),
    participant_label = factor(
      participant_label,
      levels = rev(participant_label)
    )
  )

transition_order <- transition_wide |>
  dplyr::select(participant_id, participant_label)

transition_long <- transition_obj$assignment_transitions_long |>
  dplyr::mutate(participant_id = as.character(participant_id)) |>
  dplyr::left_join(transition_order, by = "participant_id") |>
  dplyr::mutate(
    analysis = partition_factor(as.character(analysis)),
    group = group_factor(cluster)
  )

p6_heat <- ggplot2::ggplot(
  transition_long,
  ggplot2::aes(analysis, participant_label, fill = group)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.35) +
  ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
  ggplot2::labs(
    title = "Participant-level assignment changes",
    subtitle = "Rows are ordered by baseline group; sex is shown in parentheses",
    x = NULL,
    y = "Participant",
    fill = "Aligned group"
  ) +
  theme_pub(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1))

ari_plot_df <- transition_obj$ari_matrix_long |>
  dplyr::mutate(
    partition_1 = dplyr::recode(
      partition_1,
      Sex_residualized = "Sex-residualized",
      Clean_strip_residualized = "Clean-strip-residualized"
    ),
    partition_2 = dplyr::recode(
      partition_2,
      Sex_residualized = "Sex-residualized",
      Clean_strip_residualized = "Clean-strip-residualized"
    ),
    partition_1 = partition_factor(partition_1),
    partition_2 = partition_factor(partition_2)
  )

p6_ari <- ggplot2::ggplot(
  ari_plot_df,
  ggplot2::aes(partition_1, partition_2, fill = ARI)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = fmt_num(ARI, 3)), size = 3.3) +
  ggplot2::scale_fill_gradient(limits = c(0, 1)) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Pairwise partition agreement",
    x = NULL,
    y = NULL,
    fill = "ARI"
  ) +
  theme_pub(8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

cluster_size_plot <- transition_obj$cluster_size_summary |>
  dplyr::mutate(
    partition = dplyr::recode(
      partition,
      Sex_residualized = "Sex-residualized",
      Clean_strip_residualized = "Clean-strip-residualized"
    ),
    partition = partition_factor(partition),
    group = group_factor(cluster)
  )

p6_sizes <- ggplot2::ggplot(
  cluster_size_plot,
  ggplot2::aes(partition, n, fill = group)
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = 0:25) +
  ggplot2::labs(
    title = "Group sizes across the three partitions",
    x = NULL,
    y = "Participants",
    fill = "Aligned group"
  ) +
  theme_pub(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

Figure6 <- (
  p6_heat | (p6_ari / p6_sizes)
) +
  patchwork::plot_layout(widths = c(1.45, 1)) +
  patchwork::plot_annotation(tag_levels = "A")

fig6_paths <- save_figure(
  Figure6,
  "Figure6_partition_assignment_transitions",
  width = 14,
  height = 10
)


# ============================================================
# SUPPLEMENTARY FIGURES
# ============================================================


# ------------------------------------------------------------
# Figure S1: repeated-subsampling consensus matrix
# ------------------------------------------------------------

consensus_order <- names(
  sort(
    cluster_baseline * 1000 +
      suppressWarnings(as.numeric(names(cluster_baseline)))
  )
)

if (anyNA(consensus_order)) {
  consensus_order <- names(cluster_baseline)[order(cluster_baseline)]
}

consensus_df <- as.data.frame(
  as.table(robustness$consensus_matrix),
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::rename(
    participant_1 = Var1,
    participant_2 = Var2,
    consensus = Freq
  ) |>
  dplyr::mutate(
    participant_1 = factor(participant_1, levels = consensus_order),
    participant_2 = factor(participant_2, levels = rev(consensus_order))
  )

pS1 <- ggplot2::ggplot(
  consensus_df,
  ggplot2::aes(participant_1, participant_2, fill = consensus)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(limits = c(0, 1)) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Repeated 80% subsampling consensus matrix",
    x = "Participant",
    y = "Participant",
    fill = "Consensus"
  ) +
  theme_pub(8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 6.5)
  )

figS1_paths <- save_figure(
  pS1,
  "FigureS1_consensus_heatmap",
  width = 8.5,
  height = 7.5
)


# ------------------------------------------------------------
# Figure S2: cluster-specific bootstrap Jaccard distributions
# ------------------------------------------------------------

jaccard_df <- robustness$bootstrap_jaccard |>
  dplyr::mutate(group = group_factor(baseline_cluster))

pS2 <- ggplot2::ggplot(
  jaccard_df,
  ggplot2::aes(group, Jaccard, fill = group)
) +
  ggplot2::geom_violin(trim = TRUE, linewidth = 0.4) +
  ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA) +
  ggplot2::stat_summary(fun = mean, geom = "point", size = 2.3) +
  ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 1.02)) +
  ggplot2::labs(
    title = "Cluster-specific participant-bootstrap stability",
    x = NULL,
    y = "Maximum Jaccard similarity",
    fill = "Exploratory group"
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "none")

figS2_paths <- save_figure(
  pS2,
  "FigureS2_bootstrap_group_Jaccard",
  width = 7,
  height = 5.5
)


# ------------------------------------------------------------
# Figure S3: descriptive raw spectra by baseline group
# ------------------------------------------------------------

raw_group_long <- raw_long |>
  dplyr::mutate(
    group = group_factor(cluster_baseline[participant_id])
  )

raw_group_summary <- raw_group_long |>
  dplyr::group_by(group, window, wavenumber) |>
  dplyr::summarise(
    q25 = stats::quantile(absorbance, 0.25),
    median = stats::median(absorbance),
    q75 = stats::quantile(absorbance, 0.75),
    .groups = "drop"
  )

pS3 <- ggplot2::ggplot(
  raw_group_summary,
  ggplot2::aes(wavenumber, median, color = group, fill = group)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = q25, ymax = q75),
    alpha = 0.18,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::scale_x_reverse() +
  ggplot2::scale_color_manual(values = group_cols, drop = FALSE) +
  ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
  ggplot2::facet_grid(
    group ~ window,
    scales = "free",
    space = "free_x"
  ) +
  ggplot2::labs(
    title = "Descriptive raw absorbance spectra by exploratory group",
    subtitle = "Lines show medians; ribbons show interquartile ranges",
    x = expression("Wavenumber (cm"^{-1} * ")"),
    y = "Absorbance",
    color = NULL,
    fill = NULL
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "none")

figS3_paths <- save_figure(
  pS3,
  "FigureS3_descriptive_group_spectra",
  width = 11,
  height = 8.5
)


# ------------------------------------------------------------
# Figure S4: invariance across independently defined strip directions
# ------------------------------------------------------------

reference_names <- strip_obj$reference_residualization_summary$reference

pair_ari_lookup <- function(a, b) {
  if (identical(a, b)) return(1)

  z <- strip_obj$reference_residualization_pairwise_ARI |>
    dplyr::filter(
      (reference_1 == a & reference_2 == b) |
        (reference_1 == b & reference_2 == a)
    )

  if (nrow(z) != 1L) {
    stop("Could not recover a unique pairwise strip-direction ARI.")
  }

  as.numeric(z$ARI[[1]])
}

strip_ari_grid <- tidyr::expand_grid(
  reference_1 = reference_names,
  reference_2 = reference_names
) |>
  dplyr::rowwise() |>
  dplyr::mutate(ARI = pair_ari_lookup(reference_1, reference_2)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    label_1 = factor(
      clean_reference_label(reference_1),
      levels = clean_reference_label(reference_names)
    ),
    label_2 = factor(
      clean_reference_label(reference_2),
      levels = rev(clean_reference_label(reference_names))
    )
  )

pS4 <- ggplot2::ggplot(
  strip_ari_grid,
  ggplot2::aes(label_1, label_2, fill = ARI)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = fmt_num(ARI, 3)), size = 3) +
  ggplot2::scale_fill_gradient(limits = c(0, 1)) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Agreement among partitions after removal of different clean-strip directions",
    x = NULL,
    y = NULL,
    fill = "ARI"
  ) +
  theme_pub(8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

figS4_paths <- save_figure(
  pS4,
  "FigureS4_clean_strip_direction_invariance",
  width = 8.5,
  height = 7.5
)


# ============================================================
# OUTPUT MANIFEST AND SESSION INFORMATION
# ============================================================

main_table_paths <- file.path(
  tables_dir,
  c(
    "Table1_pediatric_cohort_characteristics.csv",
    "Table2_exploratory_partition_stability.csv",
    "Table3_sex_and_substrate_sensitivity.csv"
  )
)

supplement_table_paths <- file.path(
  supplement_dir,
  c(
    "TableS1_cluster_stability_detailed.csv",
    "TableS2_PC1_PC3_sex_associations_BH.csv",
    "TableS3_clean_strip_technical_reproducibility.csv",
    "TableS4_clean_strip_direction_sensitivity.csv",
    "TableS5_participant_assignment_transitions.csv",
    "TableS6_pairwise_partition_agreement.csv",
    "TableS7_PCA_eigenvector_loadings_PC1_PC3.csv"
  )
)

main_figure_paths <- c(
  fig1_paths,
  fig2_paths,
  fig3_paths,
  fig4_paths,
  fig5_paths,
  fig6_paths
)

supplement_figure_paths <- c(
  figS1_paths,
  figS2_paths,
  figS3_paths,
  figS4_paths
)

manifest <- dplyr::bind_rows(
  tibble::tibble(category = "main_table", path = main_table_paths),
  tibble::tibble(category = "supplement_table", path = supplement_table_paths),
  tibble::tibble(category = "main_figure", path = main_figure_paths),
  tibble::tibble(category = "supplement_figure", path = supplement_figure_paths)
) |>
  dplyr::mutate(
    exists = file.exists(path),
    size_bytes = ifelse(exists, file.info(path)$size, NA_real_),
    md5 = ifelse(exists, unname(tools::md5sum(path)), NA_character_)
  )

if (!all(manifest$exists)) {
  print(manifest[!manifest$exists, ], n = Inf, width = Inf)
  stop("One or more publication outputs were not created.")
}

if (any(manifest$size_bytes <= 0, na.rm = TRUE)) {
  stop("One or more publication outputs are empty.")
}

readr::write_csv(
  manifest,
  file.path(tables_dir, "07_output_manifest.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(tables_dir, "sessionInfo_tables_figures.txt")
)

cat("\n========================================\n")
cat("PEDIATRIC-ONLY TABLES AND FIGURES COMPLETE\n")
cat("========================================\n")
cat("Main tables: ", length(main_table_paths), "\n", sep = "")
cat("Main figures: ", length(main_figure_paths) / 2, " (PNG + PDF)\n", sep = "")
cat("Supplementary tables: ", length(supplement_table_paths), "\n", sep = "")
cat(
  "Supplementary figures: ",
  length(supplement_figure_paths) / 2,
  " (PNG + PDF)\n",
  sep = ""
)
cat("Manifest: tables/07_output_manifest.csv\n")
cat("========================================\n")
