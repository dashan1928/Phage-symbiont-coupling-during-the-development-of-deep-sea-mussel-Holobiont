#!/usr/bin/env Rscript
# ============================================================
# run_all_analyses.R
# Master script: source this to reproduce every figure and table.
# ============================================================

scripts <- c(
  "00_setup_and_load.R",
  "01_Fig1_vOTU_overview.R",
  "02_Fig2_diversity.R",
  "03_Fig3_taxonomy_DEseq.R",
  "04_Fig4_host_prediction.R",
  "05_Fig5_virus_host_coupling.R",
  "06_Fig6_AMG.R",
  "07_supplementary_figures.R",
  "08_supplementary_tables.R",
  "09_graphical_abstract.R"
)

for (s in scripts) {
  cat(sprintf("\n========== Running %s ==========\n", s))
  source(s, echo = FALSE)
}

cat("\nAll figures and tables generated. See ./figures and ./tables.\n")
