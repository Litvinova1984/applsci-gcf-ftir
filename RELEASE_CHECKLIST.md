# Release checklist

Before pushing the repository to GitHub:

1. Open the repository root in R/RStudio.
2. Run:

```r
install.packages("renv")  # only if renv is not installed
renv::restore()
source("R/run_all_public.R")
source("R/render_revision_report.R")
```

3. Confirm the console ends with the script-08 checks:
   - baseline cluster sizes `15/4/6`;
   - baseline 3-PC cumulative variance ≈ `89.282%`;
   - adult retained variables `793`;
   - adult mapping `PC1->PC2; PC2->PC1; PC4->PC4; PC5->PC5`;
   - domain map `D1/D2->PC2; D3/D4->PC4; D5->PC5`;
4. Confirm `analysis/revision_report.html` renders without errors.
5. Visually inspect regenerated Figures 1–7 and S1–S5; S1/S3/S4/S5 should be color figures and Figure 6 should have no internal plot title.
6. Push only after the end-to-end run succeeds.
