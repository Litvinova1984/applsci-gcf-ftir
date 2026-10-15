# ============================================================
# 01_import_preprocess.R
#
# Pediatric GCF FTIR — import and preprocessing
#
# Source:
#   data/pediatric/raw/OPUS.zip
#
# Pipeline:
#   OPUS files
#       -> opusreader2::read_opus()
#       -> AB absorbance block
#       -> select analytical spectra
#       -> retain 3400–2800 and 1800–870 cm^-1
#       -> active pediatric samples (n = 25)
#       -> X_active
#       -> row-wise SNV
#       -> X_snv
#
# No smoothing
# No baseline correction
# No derivatives
# No interpolation
# No feature-wise scaling
#
# PCA and clustering are handled in later scripts.
# ============================================================


# ------------------------------------------------------------
# 0. Project setup
# ------------------------------------------------------------

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the applsci_revision project root.")
}

source(file.path("R", "00_utils.R"))


required_pkgs <- c(
  "opusreader2",
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


cat("R version: ", R.version.string, "\n", sep = "")
cat(
  "opusreader2 version: ",
  as.character(utils::packageVersion("opusreader2")),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

zip_path <- file.path(
  "data",
  "pediatric",
  "raw",
  "OPUS.zip"
)

processed_dir <- file.path(
  "data",
  "pediatric",
  "processed"
)

if (!file.exists(zip_path)) {
  stop(
    "OPUS.zip not found at:\n",
    normalizePath(zip_path, mustWork = FALSE)
  )
}

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n========================================\n")
cat("SOURCE ARCHIVE\n")
cat("========================================\n")

cat(
  normalizePath(zip_path),
  "\n"
)

cat(
  "MD5: ",
  unname(tools::md5sum(zip_path)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Analytical source files
# ------------------------------------------------------------
#
# Selected set:
#   25 active pediatric samples
#   3 pediatric reference spectra
#   1 gingivitis reference
#   1 clean-substrate reference
#
# Alternative/unused files present in OPUS.zip
# (21_2, 22_1, 23_1, 26, Av, 17*.dpt)
# are deliberately excluded.

active_file_map <- c(
  "1"  = "1.0",
  "2"  = "2.spc.0",
  "3"  = "3.spc.0",
  "4"  = "4.spc.0",
  "5"  = "5.spc.0",
  "6"  = "6.spc.0",
  "7"  = "7.spc.0",
  "8"  = "8.spc.0",
  "9"  = "9.spc.0",
  "10" = "10.spc.0",
  "11" = "11.spc.0",
  "12" = "12.spc.0",
  "13" = "13.spc.0",
  "14" = "14.spc.0",
  "15" = "15.spc.0",
  "16" = "16.spc.0",
  "17" = "17.spc.0",
  "18" = "18.spc.0",
  "19" = "19.spc.0",
  "20" = "20.spc.0",
  "21" = "21_1.spc.0",
  "22" = "22_2.spc.0",
  "23" = "23_2.spc.0",
  "24" = "24.spc.0",
  "25" = "25.spc.0"
)

supplementary_file_map <- c(
  "etalon" = "etalon.spc.0",
  "G"      = "G.spc.0",
  "N-1"    = "N-1.spc.0",
  "N-2"    = "N-2.spc.0",
  "N"      = "N.spc.0"
)

analysis_file_map <- c(
  active_file_map,
  supplementary_file_map
)

analysis_filenames <- unname(
  analysis_file_map
)


if (length(active_file_map) != 25L) {
  stop("Expected 25 active pediatric source files.")
}

if (length(supplementary_file_map) != 5L) {
  stop("Expected 5 supplementary/reference source files.")
}

if (length(analysis_file_map) != 30L) {
  stop("Expected 30 source files in the selected analytical set.")
}


# ------------------------------------------------------------
# 3. Inspect archive contents
# ------------------------------------------------------------

zip_inventory <- utils::unzip(
  zip_path,
  list = TRUE
)

# Ignore directories and macOS resource-fork files.
zip_members <- zip_inventory$Name[
  !grepl("/$", zip_inventory$Name) &
    !grepl("^__MACOSX/", zip_inventory$Name)
]

zip_basenames <- basename(
  zip_members
)


cat("\n========================================\n")
cat("ARCHIVE CONTENTS\n")
cat("========================================\n")

cat(
  "Files in archive: ",
  length(zip_members),
  "\n",
  sep = ""
)

print(zip_basenames)


# ------------------------------------------------------------
# 4. Verify selected files are present exactly once
# ------------------------------------------------------------

missing_files <- setdiff(
  analysis_filenames,
  zip_basenames
)

if (length(missing_files) > 0) {
  stop(
    "Required source files missing from OPUS.zip:\n",
    paste(missing_files, collapse = "\n")
  )
}

duplicate_matches <- vapply(
  analysis_filenames,
  function(x) sum(zip_basenames == x),
  integer(1)
)

if (any(duplicate_matches != 1L)) {
  stop(
    "At least one required source filename does not occur exactly once."
  )
}


selected_members <- vapply(
  analysis_filenames,
  function(x) {
    zip_members[zip_basenames == x]
  },
  character(1)
)


cat(
  "\nSelected analytical set found: ",
  length(selected_members),
  " / ",
  length(analysis_filenames),
  " files\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Record excluded archive files
# ------------------------------------------------------------

excluded_members <- setdiff(
  zip_members,
  selected_members
)

cat("\nFiles present in archive but not used in the analytical set:\n")

print(
  basename(excluded_members)
)


# ------------------------------------------------------------
# 6. Extract selected files to a temporary directory
# ------------------------------------------------------------
#
# Nothing is added to or modified inside data/pediatric/raw.

tmp_dir <- tempfile(
  pattern = "pediatric_OPUS_"
)

dir.create(
  tmp_dir,
  recursive = TRUE
)

utils::unzip(
  zip_path,
  files = selected_members,
  exdir = tmp_dir,
  overwrite = TRUE
)


extracted_paths <- file.path(
  tmp_dir,
  selected_members
)

if (!all(file.exists(extracted_paths))) {
  stop("Not all selected OPUS files were extracted successfully.")
}


# ------------------------------------------------------------
# 7. OPUS reader
# ------------------------------------------------------------

read_one_opus <- function(path) {

  obj <- opusreader2::read_opus(path)

  one <- obj[[1]]

  tibble::tibble(
    source_id = tools::file_path_sans_ext(
      basename(path)
    ),
    wavenumber = one$ab$wavenumbers,
    absorbance = as.numeric(
      one$ab$data[1, ]
    )
  )
}


# ------------------------------------------------------------
# 8. Read selected spectra
# ------------------------------------------------------------

raw_outputs <- vector(
  "list",
  length(extracted_paths)
)

names(raw_outputs) <- names(
  analysis_file_map
)

read_status <- vector(
  "list",
  length(extracted_paths)
)


for (i in seq_along(extracted_paths)) {

  id <- names(analysis_file_map)[i]
  path <- extracted_paths[i]
  filename <- analysis_filenames[i]

  cat(
    "Reading ",
    filename,
    " ... ",
    sep = ""
  )

  z <- try(
    read_one_opus(path),
    silent = TRUE
  )

  if (inherits(z, "try-error")) {

    cat("FAILED\n")

    read_status[[i]] <- tibble::tibble(
      sample_id = id,
      source_file = filename,
      readable = FALSE,
      error = as.character(z)[1]
    )

  } else {

    cat(
      "OK (",
      nrow(z),
      " points)\n",
      sep = ""
    )

    raw_outputs[[i]] <- z

    read_status[[i]] <- tibble::tibble(
      sample_id = id,
      source_file = filename,
      readable = TRUE,
      error = NA_character_
    )
  }
}


read_status <- dplyr::bind_rows(
  read_status
)


cat("\n========================================\n")
cat("OPUS IMPORT CHECK\n")
cat("========================================\n")

print(
  read_status,
  n = Inf
)


if (any(!read_status$readable)) {

  stop(
    paste0(
      "\nAt least one selected OPUS file could not be read.\n",
      "Check the source archive and the opusreader2/R environment."
    )
  )
}


# ------------------------------------------------------------
# 9. Add canonical sample IDs and roles
# ------------------------------------------------------------

spectra_list <- vector(
  "list",
  length(raw_outputs)
)

for (i in seq_along(raw_outputs)) {

  z <- raw_outputs[[i]]

  id <- names(analysis_file_map)[i]

  role_i <- if (id %in% names(active_file_map)) {
    "Active"
  } else {
    "Supplementary"
  }

  z$sample_id <- id
  z$source_file <- analysis_filenames[i]
  z$role <- role_i

  spectra_list[[i]] <- z
}

spectra_long <- dplyr::bind_rows(
  spectra_list
)

unlink(
  tmp_dir,
  recursive = TRUE,
  force = TRUE
)


# ------------------------------------------------------------
# 10. Basic spectral integrity checks
# ------------------------------------------------------------

if (anyNA(spectra_long$wavenumber)) {
  stop("NA wavenumbers detected.")
}

if (anyNA(spectra_long$absorbance)) {
  stop("NA absorbance values detected.")
}

if (any(!is.finite(spectra_long$wavenumber))) {
  stop("Non-finite wavenumbers detected.")
}

if (any(!is.finite(spectra_long$absorbance))) {
  stop("Non-finite absorbance values detected.")
}


points_per_spectrum <- spectra_long |>
  dplyr::count(
    sample_id,
    name = "n_points"
  )

if (dplyr::n_distinct(points_per_spectrum$n_points) != 1L) {
  stop("Source spectra do not contain the same number of spectral points.")
}


# ------------------------------------------------------------
# 11. Active pediatric cohort
# ------------------------------------------------------------

spectra_active <- spectra_long |>
  dplyr::filter(
    role == "Active"
  )


active_ids_found <- unique(
  spectra_active$sample_id
)

if (!setequal(
  active_ids_found,
  as.character(1:25)
)) {
  stop(
    "Active cohort is not exactly sample IDs 1:25."
  )
}


# ------------------------------------------------------------
# 12. Retain analytical windows
# ------------------------------------------------------------
#
# Analytical windows:
#   3400–2800 cm^-1
#   1800–870 cm^-1
#
# The retained instrument grid is approximately
# 3399.01 ... 2800.xx and 1799.xx ... 870.50 cm^-1.
#
# No interpolation is performed.

spectra_active_window <- spectra_active |>
  dplyr::filter(
    retain_analysis_windows(
      wavenumber
    )
  )


# ------------------------------------------------------------
# 13. Deterministic sample order
# ------------------------------------------------------------
#
# Samples are ordered lexicographically by source filename.

active_manifest <- tibble::tibble(
  sample_id = names(active_file_map),
  source_file = unname(active_file_map)
) |>
  dplyr::arrange(
    source_file
  )

active_order <- active_manifest$sample_id


cat("\n========================================\n")
cat("ACTIVE MATRIX ROW ORDER\n")
cat("========================================\n")

print(active_order)


# ------------------------------------------------------------
# 14. Extract one retained active spectrum
# ------------------------------------------------------------

get_active_spectrum <- function(id) {

  z <- spectra_active_window |>
    dplyr::filter(
      sample_id == id
    ) |>
    dplyr::arrange(
      dplyr::desc(wavenumber)
    )

  if (nrow(z) == 0L) {
    stop(
      "No retained spectral data for active sample ",
      id
    )
  }

  z
}


# ------------------------------------------------------------
# 15. Establish common wavenumber grid
# ------------------------------------------------------------

first_spectrum <- get_active_spectrum(
  active_order[1]
)

wn_ref <- as.numeric(
  first_spectrum$wavenumber
)


for (id in active_order) {

  z <- get_active_spectrum(id)

  wn_i <- as.numeric(
    z$wavenumber
  )

  if (length(wn_i) != length(wn_ref)) {

    stop(
      "Different retained grid length for sample ",
      id,
      ": ",
      length(wn_i),
      " vs ",
      length(wn_ref)
    )
  }

  if (!isTRUE(
    all.equal(
      wn_i,
      wn_ref,
      tolerance = 1e-12
    )
  )) {

    stop(
      "Wavenumber-grid mismatch for sample ",
      id,
      ". No interpolation was performed."
    )
  }
}


# ------------------------------------------------------------
# 16. Construct X_active
# ------------------------------------------------------------

X_active <- do.call(
  rbind,
  lapply(
    active_order,
    function(id) {

      as.numeric(
        get_active_spectrum(id)$absorbance
      )
    }
  )
)

rownames(X_active) <- active_order

colnames(X_active) <- make_wn_names(
  wn_ref
)


# ------------------------------------------------------------
# 17. Structural validation
# ------------------------------------------------------------

if (!identical(
  dim(X_active),
  c(25L, 793L)
)) {

  stop(
    "Expected X_active dimensions 25 x 793; obtained ",
    paste(dim(X_active), collapse = " x ")
  )
}

if (any(!is.finite(X_active))) {
  stop("Non-finite values in X_active.")
}

assert_no_duplicates(
  rownames(X_active),
  "sample IDs"
)

assert_no_duplicates(
  colnames(X_active),
  "spectral column names"
)


cat("\n========================================\n")
cat("X_ACTIVE\n")
cat("========================================\n")

cat(
  "Dimensions: ",
  paste(dim(X_active), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "First wavenumbers: ",
  paste(head(wn_ref), collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Last wavenumbers: ",
  paste(tail(wn_ref), collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Mean: ",
  mean(X_active),
  "\n",
  sep = ""
)

cat(
  "SD: ",
  stats::sd(as.numeric(X_active)),
  "\n",
  sep = ""
)

cat(
  "Range: ",
  paste(range(X_active), collapse = " to "),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Row-wise SNV
# ------------------------------------------------------------
#
# Each participant spectrum is normalized independently:
#   (x - mean(x)) / sd(x)

X_snv <- snv_matrix(
  X_active
)


snv_row_mean <- rowMeans(
  X_snv
)

snv_row_sd <- apply(
  X_snv,
  1,
  stats::sd
)


if (max(abs(snv_row_mean)) > 1e-10) {
  stop("SNV check failed: row means are not approximately zero.")
}

if (max(abs(snv_row_sd - 1)) > 1e-10) {
  stop("SNV check failed: row SDs are not approximately one.")
}


# ------------------------------------------------------------
# 19. Save processed analytical objects
# ------------------------------------------------------------

saveRDS(
  X_active,
  file.path(
    processed_dir,
    "pediatric_X_active_absorbance.rds"
  )
)

saveRDS(
  X_snv,
  file.path(
    processed_dir,
    "pediatric_X_snv.rds"
  )
)

saveRDS(
  wn_ref,
  file.path(
    processed_dir,
    "pediatric_wavenumbers.rds"
  )
)

saveRDS(
  spectra_long,
  file.path(
    processed_dir,
    "pediatric_source_spectra_long.rds"
  )
)


# ------------------------------------------------------------
# 20. Portable CSV outputs
# ------------------------------------------------------------

readr::write_csv(
  data.frame(
    sample_id = rownames(X_active),
    X_active,
    check.names = FALSE
  ),
  file.path(
    processed_dir,
    "pediatric_X_active_absorbance.csv"
  )
)

readr::write_csv(
  data.frame(
    sample_id = rownames(X_snv),
    X_snv,
    check.names = FALSE
  ),
  file.path(
    processed_dir,
    "pediatric_X_snv.csv"
  )
)


# ------------------------------------------------------------
# 21. Source manifest
# ------------------------------------------------------------

source_manifest <- tibble::tibble(
  sample_id = names(analysis_file_map),
  source_file = unname(analysis_file_map),
  role = c(
    rep("Active", length(active_file_map)),
    rep(
      "Supplementary",
      length(supplementary_file_map)
    )
  )
)

readr::write_csv(
  source_manifest,
  file.path(
    processed_dir,
    "pediatric_source_manifest.csv"
  )
)


readr::write_csv(
  read_status,
  file.path(
    processed_dir,
    "pediatric_OPUS_import_status.csv"
  )
)


# ------------------------------------------------------------
# 22. Archive inventory
# ------------------------------------------------------------

archive_inventory <- tibble::tibble(
  archive_member = zip_members,
  source_file = basename(zip_members),
  used_in_analytical_set =
    zip_members %in% selected_members
)

readr::write_csv(
  archive_inventory,
  file.path(
    processed_dir,
    "OPUS_archive_inventory.csv"
  )
)


# ------------------------------------------------------------
# 23. Numerical preprocessing audit
# ------------------------------------------------------------

preprocessing_audit <- tibble::tibble(

  item = c(
    "Selected source spectra",
    "Active pediatric spectra",
    "Supplementary/reference spectra",
    "Retained spectral variables",
    "Highest retained wavenumber",
    "Lowest retained wavenumber",
    "Mean X_active",
    "SD X_active",
    "Maximum absolute SNV row mean",
    "Maximum absolute deviation of SNV row SD from 1"
  ),

  value = c(
    length(analysis_file_map),
    length(active_file_map),
    length(supplementary_file_map),
    ncol(X_active),
    max(wn_ref),
    min(wn_ref),
    mean(X_active),
    stats::sd(as.numeric(X_active)),
    max(abs(snv_row_mean)),
    max(abs(snv_row_sd - 1))
  )
)

readr::write_csv(
  preprocessing_audit,
  file.path(
    processed_dir,
    "preprocessing_audit.csv"
  )
)


# ------------------------------------------------------------
# 24. Reproducibility environment
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    processed_dir,
    "sessionInfo_import_preprocess.txt"
  )
)


# ------------------------------------------------------------
# 25. Final report
# ------------------------------------------------------------

cat("\n\n")
cat("========================================\n")
cat("PEDIATRIC PREPROCESSING COMPLETE\n")
cat("========================================\n")

cat(
  "Selected source files:    ",
  length(analysis_file_map),
  "\n",
  sep = ""
)

cat(
  "Active spectra:           ",
  nrow(X_active),
  "\n",
  sep = ""
)

cat(
  "Spectral variables:       ",
  ncol(X_active),
  "\n",
  sep = ""
)

cat(
  "Highest wavenumber:       ",
  max(wn_ref),
  " cm^-1\n",
  sep = ""
)

cat(
  "Lowest wavenumber:        ",
  min(wn_ref),
  " cm^-1\n",
  sep = ""
)

cat(
  "X_active mean:            ",
  mean(X_active),
  "\n",
  sep = ""
)

cat(
  "X_active SD:              ",
  stats::sd(as.numeric(X_active)),
  "\n",
  sep = ""
)

cat(
  "Max |SNV row mean|:       ",
  max(abs(snv_row_mean)),
  "\n",
  sep = ""
)

cat(
  "Max |SNV row SD - 1|:     ",
  max(abs(snv_row_sd - 1)),
  "\n",
  sep = ""
)

cat(
  "\nOutputs written to:\n",
  normalizePath(processed_dir),
  "\n",
  sep = ""
)

cat("========================================\n")
