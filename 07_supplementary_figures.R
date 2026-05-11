# ============================================================
# 07_supplementary_figures.R
# Figures S1-S5
# ============================================================

# ── FigS1: NMDS validation ──────────────────────────────────
set.seed(42)
nmds_res <- vegan::metaMDS(comm, distance = "bray", k = 2,
                           trymax = 100, trace = FALSE)

nmds_df <- as.data.frame(nmds_res$points) %>%
  setNames(c("NMDS1","NMDS2")) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels))

figS1 <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = Stage, shape = Stage)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(type = "t", linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = stage_pal) +
  scale_shape_manual(values = c(Early = 16, Middle = 17, Late = 15)) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Stress = %.3f", nmds_res$stress),
           hjust = 1.05, vjust = 1.1, size = 3.5,
           fill = "white", color = "grey30", label.size = 0.3) +
  labs(title = "Figure S1: NMDS Validation",
       subtitle = "Bray-Curtis dissimilarity (k=2)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey45"))

ggsave(file.path(FIG_DIR, "FigS1_NMDS.pdf"), figS1,
       width = 7, height = 6, device = cairo_pdf)
cat("FigS1 saved\n")

# ── FigS2: CheckV quality ───────────────────────────────────
checkv_all <- read_tsv(file.path(DATA_DIR, "merged_checkv_quality.tsv"),
                       show_col_types = FALSE) %>%
  dplyr::rename(SampleID = Sample) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  filter(!is.na(Stage)) %>%
  mutate(Stage = factor(Stage, levels = stage_levels),
         quality = factor(checkv_quality,
                          levels = c("High-quality","Medium-quality",
                                     "Low-quality","Not-determined")))

per_sample_n <- checkv_all %>% dplyr::count(SampleID, Stage)
kw_yield <- kruskal.test(n ~ Stage, data = per_sample_n)

figS2a <- ggplot(per_sample_n, aes(Stage, n, fill = Stage, color = Stage)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = 0.55, width = 0.52) +
  geom_jitter(width = 0.18, size = 2, alpha = 0.75,
              shape = 21, stroke = 0.25, color = "white") +
  annotate("text", x = 2, y = Inf,
           label = sprintf("KW p = %.3f", kw_yield$p.value),
           vjust = -0.4, size = 3.4, fontface = "italic", color = "grey40") +
  scale_fill_manual(values = stage_pal, guide = "none") +
  scale_color_manual(values = stage_pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "vOTU count per sample",
       title = "Per-Sample vOTU Yield") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank())

qual_cols <- c("High-quality" = "#1B5E20", "Medium-quality" = "#66BB6A",
               "Low-quality" = "#FFA726", "Not-determined" = "#CFD8DC")
qual_cnt <- checkv_all %>% dplyr::count(quality) %>%
  mutate(pct = n / sum(n) * 100,
         quality = factor(quality, levels = names(qual_cols)))

figS2b <- ggplot(qual_cnt, aes(x = 1, y = n, fill = quality)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%s\n%d (%.0f%%)", quality, n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.0, color = "white", fontface = "bold", lineheight = 0.85) +
  coord_polar("y") + xlim(0, 1.5) +
  scale_fill_manual(values = qual_cols, guide = "none") +
  labs(title = "CheckV Quality Distribution") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

figS2 <- patchwork::wrap_plots(figS2a, figS2b, nrow = 1, widths = c(1.2, 1)) +
  patchwork::plot_annotation(
    title = "Figure S2: Viral Contig Quality Assessment",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(FIG_DIR, "FigS2_quality.pdf"), figS2,
       width = 11, height = 5, device = cairo_pdf)
cat("FigS2 saved\n")

# ── FigS3: bacterial MAG composition ────────────────────────
mag_stage_comp <- mag_tpm_long %>%
  mutate(phylum_lab = case_when(
    phylum %in% c("Pseudomonadota","Campylobacterota","Actinomycetota",
                  "Bacillota","Spirochaetota") ~ phylum,
    TRUE ~ "Other")) %>%
  group_by(SampleID, Stage, phylum_lab) %>%
  summarise(tot = sum(MAG_TPM, na.rm = TRUE), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(prop = tot / sum(tot) * 100,
         Stage = factor(Stage, levels = stage_levels)) %>% ungroup()

phy_cols <- c(Pseudomonadota = "#1565C0", Campylobacterota = "#C62828",
              Actinomycetota = "#2E7D32", Bacillota = "#F9A825",
              Spirochaetota = "#6A1B9A", Other = "#CFD8DC")

figS3 <- ggplot(
  mag_stage_comp %>%
    mutate(SampleID = factor(SampleID, levels = sample_order_vec),
           phylum_lab = factor(phylum_lab, levels = names(phy_cols))),
  aes(SampleID, prop, fill = phylum_lab)) +
  geom_col(width = 0.92) +
  facet_grid(~ Stage, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phy_cols, name = "Bacterial Phylum") +
  labs(x = NULL, y = "Relative Abundance (% TPM)",
       title    = "Figure S3: Bacterial MAG Community Composition",
       subtitle = "Phylum-level | 23 MAGs across 48 samples") +
  theme_bw(base_size = 10) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 6.5),
        strip.text       = element_text(face = "bold", size = 11, color = "white"),
        strip.background = element_rect(fill = "#37474F"),
        plot.title       = element_text(face = "bold"),
        panel.grid.major.x = element_blank())

ggsave(file.path(FIG_DIR, "FigS3_MAG_composition.pdf"), figS3,
       width = 14, height = 5, device = cairo_pdf)
cat("FigS3 saved\n")

# ── FigS4: MA plots ─────────────────────────────────────────
ma_df <- deseq_all %>%
  mutate(comparison = factor(comparison,
                             levels = c("Early vs Middle","Early vs Late","Middle vs Late")),
         sig = case_when(
           padj < 0.05 & log2FoldChange >  1 ~ "Up in B",
           padj < 0.05 & log2FoldChange < -1 ~ "Up in A",
           TRUE ~ "NS"))

figS4 <- ggplot(ma_df, aes(log10(baseMean + 1), log2FoldChange,
                           color = sig, size = sig)) +
  geom_point(alpha = 0.55) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_color_manual(values = c("Up in A" = "#0072B2", "Up in B" = "#D55E00", NS = "grey78"),
                     name = NULL) +
  scale_size_manual(values = c("Up in A" = 1.8, "Up in B" = 1.8, NS = 0.8),
                    guide = "none") +
  labs(x = "log10(mean normalised count + 1)", y = "log2 Fold Change",
       title    = "Figure S4: MA Plots \u2014 DESeq2 Differential Abundance",
       subtitle = "Dashed lines at |LFC| = 1; BH-adjusted p < 0.05") +
  theme_bw(base_size = 11) +
  theme(plot.title       = element_text(face = "bold"),
        strip.text       = element_text(face = "bold", color = "white"),
        strip.background = element_rect(fill = "#37474F"),
        legend.position  = "top")

ggsave(file.path(FIG_DIR, "FigS4_MAplots.pdf"), figS4,
       width = 14, height = 5, device = cairo_pdf)
cat("FigS4 saved\n")

# ── FigS5: top 30 KO ────────────────────────────────────────
top30_ko <- amg %>%
  filter(rank == "B", module != "Other") %>%
  dplyr::count(ko_id, ko_name, module, sort = TRUE) %>%
  slice_head(n = 30) %>%
  mutate(ko_short = str_extract(ko_name, "^[^;\\[]+") %>%
         str_trim() %>% str_sub(1, 38),
         label = sprintf("%s (%s)", ko_short, ko_id),
         label = factor(label, levels = rev(unique(label))),
         module = factor(module, levels = names(mod_pal)))

figS5 <- ggplot(top30_ko, aes(label, n, fill = module)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), hjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = mod_pal, name = "Module") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(x = NULL, y = "Copy number (rank-B AMGs)",
       title    = "Figure S5: Top 30 AMG KOs (DRAM-v rank B)") +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

ggsave(file.path(FIG_DIR, "FigS5_top30KO.pdf"), figS5,
       width = 11, height = 9, device = cairo_pdf)
cat("FigS5 saved\n")
