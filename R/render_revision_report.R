# ============================================================
# render_revision_report.R
#
# Render analysis/revision_report.qmd to HTML.
# Run from the project root.
# ============================================================

qmd <- file.path(
  "analysis",
  "revision_report.qmd"
)

if (!file.exists(qmd)) {
  stop("Missing ", qmd, ".")
}

quarto_bin <- Sys.which("quarto")

if (!nzchar(quarto_bin)) {
  stop(
    "Quarto executable was not found on PATH. ",
    "Open analysis/revision_report.qmd in RStudio/Posit and click Render, ",
    "or install Quarto and rerun."
  )
}

status <- system2(
  quarto_bin,
  c(
    "render",
    qmd,
    "--to",
    "html"
  )
)

if (!identical(status, 0L)) {
  stop(
    "Quarto render failed with exit status ",
    status,
    "."
  )
}

html <- file.path(
  "analysis",
  "revision_report.html"
)

if (!file.exists(html)) {
  stop(
    "Render completed but ",
    html,
    " was not found."
  )
}

cat(
  "\n========================================\n",
  "REVISION REPORT RENDERED\n",
  "========================================\n",
  "HTML: ",
  normalizePath(
    html,
    mustWork = TRUE
  ),
  "\n========================================\n",
  sep = ""
)
