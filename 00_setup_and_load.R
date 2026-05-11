# ============================================================
# 00_setup_and_load.R
# Load all data, packages, define palettes and helpers used by
# every subsequent figure script.
#
# This script must be sourced FIRST. All Fig*.R scripts assume
# the objects it creates (tpm_mat, meta_sub, host_lookup, etc.)
# already exist in the global environment.
# ============================================================

# ── Required packages ───────────────────────────────────────
required_pkgs <- c(
  "tidyverse", "readxl", "writexl",
  "vegan", "ape",
  "DESeq2", "lme4", "broom.mixed",
  "dunn.test", "ggsignif", "ggrepel",
  "patchwork", "pheatmap", "RColorBrewer", "scales"
)

missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  stop(sprintf("Install required packages first: %s",
               paste(missing, collapse = ", ")))
}
suppressPackageStartupMessages({
  invisible(lapply(required_pkgs, library, character.only = TRUE))
})

# ── Project paths (edit if needed) ──────────────────────────
DATA_DIR    <- "data"        # directory with the input tables
FIG_DIR     <- "figures"
TABLE_DIR   <- "tables"
dir.create(FIG_DIR,   showWarnings = FALSE)
dir.create(TABLE_DIR, showWarnings = FALSE)

# ── Stage palette (Okabe-Ito; colourblind-safe) ─────────────
stage_pal    <- c(Early = "#0072B2", Middle = "#009E73", Late = "#D55E00")
stage_levels <- c("Early", "Middle", "Late")

# ── AMG functional-module palette ───────────────────────────
mod_pal <- c(
  "DNA Replication/Repair"      = "#1565C0",
  "Signaling/Energy Metabolism" = "#E65100",
  "DNA Methylation"             = "#6A1B9A",
  "Biosynthesis"                = "#2E7D32",
  "Toxin-Antitoxin"             = "#C62828",
  "Electron Transport Chain"    = "#00838F",
  "Chaperones"                  = "#4E342E",
  "Transporters"                = "#F9A825",
  "Lipid Metabolism"            = "#558B2F",
  "Translation"                 = "#AD1457",
  "CRISPR-Cas"                  = "#37474F",
  "Sulfur Metabolism"           = "#FF8F00",
  "Methyl Cycle"                = "#7B1FA2",
  "Other"                       = "#90A4AE"
)

# ── Module label translation (CN -> EN) ─────────────────────
mod_map <- c(
  "其他"          = "Other",
  "DNA 复制/修复" = "DNA Replication/Repair",
  "信号/能量代谢" = "Signaling/Energy Metabolism",
  "DNA 甲基化"    = "DNA Methylation",
  "合成代谢"      = "Biosynthesis",
  "毒素-抗毒素"   = "Toxin-Antitoxin",
  "电子传递链"    = "Electron Transport Chain",
  "分子伴侣"      = "Chaperones",
  "转运蛋白"      = "Transporters",
  "脂质代谢"      = "Lipid Metabolism",
  "翻译系统"      = "Translation",
  "CRISPR-Cas"    = "CRISPR-Cas"
)

# ── Helper: significance label ──────────────────────────────
sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

# ============================================================
# Load metadata
# ============================================================
meta_sub <- read_csv(file.path(DATA_DIR, "sample_metadata.csv"),
                     show_col_types = FALSE) %>%
  mutate(Stage = factor(Stage, levels = stage_levels))
cat(sprintf("Metadata loaded: %d samples\n", nrow(meta_sub)))

# ============================================================
# Load TPM matrices
# ============================================================
tpm_raw <- read_csv(file.path(DATA_DIR, "virus_only_tpm_matrix.csv"),
                    show_col_types = FALSE)
if (!"Contig_norm" %in% colnames(tpm_raw)) {
  colnames(tpm_raw)[1] <- "Contig_norm"
}

tpm_mat <- tpm_raw %>%
  column_to_rownames("Contig_norm") %>%
  as.matrix()
tpm_mat <- tpm_mat[, meta_sub$SampleID]
cat(sprintf("vOTU TPM matrix: %d vOTUs x %d samples\n",
            nrow(tpm_mat), ncol(tpm_mat)))

mag_tpm_raw <- read_tsv(file.path(DATA_DIR, "bacteria_host_MAG_TPM.tsv"),
                        show_col_types = FALSE)
mag_tpm_mat <- mag_tpm_raw %>%
  column_to_rownames("MAG_ID") %>%
  as.matrix()
mag_tpm_mat <- mag_tpm_mat[, meta_sub$SampleID]
cat(sprintf("MAG TPM matrix:  %d MAGs x %d samples\n",
            nrow(mag_tpm_mat), ncol(mag_tpm_mat)))

# ============================================================
# Load geNomad taxonomy (for vOTUs)
# ============================================================
genomad_tax <- read_tsv(file.path(DATA_DIR, "merged_genomad_summary.tsv"),
                        show_col_types = FALSE) %>%
  dplyr::transmute(
    Contig_norm  = str_remove(seq_name, "\\|\\|.*$"),
    family_raw   = str_split_fixed(taxonomy, ";", 7)[, 7],
    class_raw    = str_split_fixed(taxonomy, ";", 7)[, 5],
    family_label = case_when(
      trimws(family_raw) != ""               ~ trimws(family_raw),
      trimws(class_raw)  == "Caudoviricetes"  ~ "Unclassified Caudoviricetes",
      trimws(class_raw)  != ""                ~ paste0("Unclassified ", trimws(class_raw)),
      TRUE                                    ~ "Unclassified"
    )
  ) %>%
  dplyr::transmute(Contig_norm, family_label) %>%
  distinct(Contig_norm, .keep_all = TRUE)
cat(sprintf("geNomad taxonomy: %d unique vOTUs annotated\n", nrow(genomad_tax)))

# ============================================================
# Load iPHoP host predictions, parse genus from full GTDB string
# ============================================================
host_pred <- read_csv(file.path(DATA_DIR, "Host_prediction_to_genus_m90.csv"),
                      show_col_types = FALSE) %>%
  dplyr::rename(
    contig_id  = Virus,
    host_genus = `Host genus`,
    confidence = `Confidence score`
  ) %>%
  mutate(contig_id = str_remove(contig_id, "\\|\\|.*$"))

# Convert GTDB-style strings to readable genus labels
host_lookup <- host_pred %>%
  dplyr::transmute(
    Contig_norm = str_remove(contig_id, "\\|\\|.*$"),
    genus_raw   = str_extract(host_genus, "(?<=g__)[^;]+") %>% str_trim(),
    order_raw   = str_extract(host_genus, "(?<=o__)[^;]+") %>% str_trim(),
    host_genus  = case_when(
      !is.na(genus_raw) &
        !str_detect(genus_raw, "^[A-Z]{2,}[0-9]+$") ~ genus_raw,
      !is.na(order_raw)                              ~
        paste0(str_remove(order_raw, "ales$"), ".(novel)"),
      TRUE                                           ~ "Unknown"
    )
  ) %>%
  dplyr::select(Contig_norm, host_genus) %>%
  distinct(Contig_norm, .keep_all = TRUE)
cat(sprintf("iPHoP host lookup: %d vOTU-host pairs\n", nrow(host_lookup)))

# ============================================================
# Load CheckV / vOTU annotation
# ============================================================
votu_annot <- read_tsv(file.path(DATA_DIR, "quality_summary.tsv"),
                       show_col_types = FALSE) %>%
  mutate(contig_id = str_remove(contig_id, "\\|\\|.*$")) %>%
  left_join(genomad_tax, by = c("contig_id" = "Contig_norm"))

# ============================================================
# Load MAG taxonomy from GTDB-Tk
# ============================================================
mag_taxonomy <- read_tsv(file.path(DATA_DIR, "gtdbtk_bac120_summary.tsv"),
                         show_col_types = FALSE) %>%
  dplyr::transmute(
    MAG_ID = user_genome,
    phylum = str_extract(classification, "(?<=p__)[^;]+") %>% str_trim(),
    order  = str_extract(classification, "(?<=o__)[^;]+") %>% str_trim(),
    family = str_extract(classification, "(?<=f__)[^;]+") %>% str_trim(),
    genus  = str_extract(classification, "(?<=g__)[^;]+") %>% str_trim(),
    genus_label = case_when(
      !is.na(genus) &
        !str_detect(genus, "^[A-Z]{2,}[0-9]+$")  ~ genus,
      !is.na(family) &
        !str_detect(family, "^[A-Z]{2,}[0-9]+$") ~ paste0(family, " (family)"),
      !is.na(order)                              ~
        paste0(str_remove(order, "ales$"), ".(novel)"),
      TRUE                                       ~ "Unknown"
    )
  )
cat(sprintf("GTDB-Tk taxonomy: %d MAGs\n", nrow(mag_taxonomy)))

# ============================================================
# Load and merge AMG annotations (DRAM-v + VIBRANT supplement)
# ============================================================
amg_dram <- read_xlsx(file.path(DATA_DIR, "dram_v_AMG_results.xlsx")) %>%
  mutate(
    scaffold = str_remove(vOTU, "^Global_Dereplicated_vOTUs\\.unique_"),
    ko_name  = kegg_hit,
    module   = mod_map[functional_category] %>% replace_na("Other"),
    source   = "DRAM-v"
  ) %>%
  dplyr::select(scaffold, vOTU, ko_id, ko_name, module,
                functional_category, host_symbiont,
                confidence, rank, source)

amg_vibrant_sup <- read_tsv(
  file.path(DATA_DIR, "VIBRANT_AMG_individuals_Global_Dereplicated_vOTUs_clean.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::rename(ko_id = `AMG KO`, ko_name = `AMG KO name`) %>%
  filter(!ko_id %in% unique(amg_dram$ko_id)) %>%
  mutate(
    scaffold = str_remove(scaffold, "\\|\\|full$"),
    vOTU     = paste0("Global_Dereplicated_vOTUs.unique_", scaffold),
    module   = case_when(
      ko_id %in% c("K00390","K01738","K21140","K00383","K04487") ~ "Sulfur Metabolism",
      ko_id == "K17398"                                          ~ "DNA Methylation",
      ko_id == "K01243"                                          ~ "Methyl Cycle",
      TRUE                                                       ~ "Other"
    ),
    functional_category = module,
    host_symbiont = NA_character_,
    confidence    = NA_character_,
    rank          = NA_character_,
    source        = "VIBRANT"
  ) %>%
  dplyr::select(scaffold, vOTU, ko_id, ko_name, module,
                functional_category, host_symbiont,
                confidence, rank, source)

amg    <- bind_rows(amg_dram, amg_vibrant_sup)
amg_hq <- amg %>% filter(rank == "B" | source == "VIBRANT")
cat(sprintf("AMG annotations: total=%d, high-confidence=%d, unique KOs=%d\n",
            nrow(amg), nrow(amg_hq), n_distinct(amg$ko_id)))

# ============================================================
# Sample ordering for heatmaps (Early -> Middle -> Late)
# ============================================================
sample_order_vec <- meta_sub %>%
  arrange(Stage) %>%
  pull(SampleID)

# Annotation data frame for heatmap columns
col_ann <- meta_sub %>%
  dplyr::select(SampleID, Stage) %>%
  column_to_rownames("SampleID")
col_ann <- col_ann[sample_order_vec, , drop = FALSE]

cat("\n=== Setup complete ===\n")
