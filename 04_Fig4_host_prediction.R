# ============================================================
# 04_Fig4_host_prediction.R
# Figure 4: host prediction donut + phage abundance + coupling preview
# ============================================================

# ── MAG TPM long table ──────────────────────────────────────
mag_tpm_long <- mag_tpm_raw %>%
  pivot_longer(-MAG_ID, names_to = "SampleID", values_to = "MAG_TPM") %>%
  filter(SampleID %in% meta_sub$SampleID) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  left_join(mag_taxonomy %>% dplyr::select(MAG_ID, genus_label, phylum),
            by = "MAG_ID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels))

# ── Panel A: host donut ─────────────────────────────────────
host_cnt_full <- host_lookup %>%
  dplyr::count(host_genus, sort = TRUE) %>%
  mutate(genus_lab = case_when(
    row_number() <= 7 & host_genus != "Unknown" ~ host_genus,
    TRUE ~ "Other")) %>%
  group_by(genus_lab) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(desc(n)) %>%
  mutate(pct = n / sum(n) * 100,
         genus_lab = factor(genus_lab,
                            levels = c(setdiff(genus_lab[order(-n)], "Other"), "Other")))

donut_cols <- setNames(
  c("#C62828", "#0072B2", brewer.pal(8, "Set2")[3:7], "grey72"),
  levels(host_cnt_full$genus_lab))

fig4a <- ggplot(host_cnt_full, aes(x = 2, y = n, fill = genus_lab)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct > 4,
                               sprintf("%s\n%.1f%%", genus_lab, pct), "")),
            position = position_stack(vjust = 0.5),
            size = 3.0, color = "white", fontface = "bold", lineheight = 0.85) +
  coord_polar("y") + xlim(0.5, 2.5) +
  scale_fill_manual(values = donut_cols, name = "Host Genus") +
  labs(title = "iPHoP Host Prediction",
       subtitle = sprintf("n = %d vOTUs | \u226590%% confidence",
                          nrow(host_lookup))) +
  theme_void() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, color = "grey45", size = 9))

# ── Panel B: MOX & CAM phage abundance per stage ────────────
mox_ids <- host_lookup %>% filter(host_genus == "Methyloprofundus") %>%
  pull(Contig_norm) %>% intersect(rownames(tpm_mat))
cam_ids <- host_lookup %>% filter(host_genus == "Campylobacter.(novel)") %>%
  pull(Contig_norm) %>% intersect(rownames(tpm_mat))

phage_tpm_df <- bind_rows(
  tibble(SampleID = colnames(tpm_mat),
         TPM_sum  = colSums(tpm_mat[mox_ids, , drop = FALSE]),
         PhageType = "Methyloprofundus phage"),
  tibble(SampleID = colnames(tpm_mat),
         TPM_sum  = colSums(tpm_mat[cam_ids, , drop = FALSE]),
         PhageType = "Campylobacter.(novel) phage")) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels),
         PhageType = factor(PhageType,
                            levels = c("Methyloprofundus phage",
                                       "Campylobacter.(novel) phage")))

kw_mox <- kruskal.test(TPM_sum ~ Stage,
                       data = filter(phage_tpm_df, PhageType == "Methyloprofundus phage"))
kw_cam <- kruskal.test(TPM_sum ~ Stage,
                       data = filter(phage_tpm_df, PhageType == "Campylobacter.(novel) phage"))

dt_mox <- dunn.test::dunn.test(
  filter(phage_tpm_df, PhageType == "Methyloprofundus phage")$TPM_sum,
  filter(phage_tpm_df, PhageType == "Methyloprofundus phage")$Stage,
  method = "BH", list = FALSE, kw = FALSE, label = FALSE)
dt_cam <- dunn.test::dunn.test(
  filter(phage_tpm_df, PhageType == "Campylobacter.(novel) phage")$TPM_sum,
  filter(phage_tpm_df, PhageType == "Campylobacter.(novel) phage")$Stage,
  method = "BH", list = FALSE, kw = FALSE, label = FALSE)

parse_comps_fn <- function(dt, alpha = 0.05) {
  keep <- dt$P.adjusted < alpha
  list(comparisons = lapply(dt$comparisons[keep],
                            function(s) trimws(strsplit(s, " - ")[[1]])),
       annotations = sig_label(dt$P.adjusted[keep]))
}

make_phage_box <- function(df, dt_obj, kw_p, title_str, fill_col) {
  pc <- parse_comps_fn(dt_obj)
  p <- ggplot(df, aes(Stage, TPM_sum)) +
    geom_boxplot(fill = fill_col, alpha = 0.35, outlier.shape = NA,
                 linewidth = 0.55, width = 0.52) +
    geom_jitter(width = 0.18, size = 2.0, alpha = 0.75,
                color = fill_col, shape = 21, fill = fill_col, stroke = 0.25) +
    scale_y_continuous(labels = scales::label_comma(),
                       expand = expansion(mult = c(0.05, 0.20))) +
    annotate("text", x = 2, y = Inf,
             label = sprintf("KW p = %.3f", kw_p),
             vjust = -0.4, size = 3.4, fontface = "italic",
             color = ifelse(kw_p < 0.05, "grey20", "grey55")) +
    labs(x = NULL, y = "Total TPM", title = title_str) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10.5),
          panel.grid.major.x = element_blank())
  if (length(pc$comparisons) > 0)
    p <- p + ggsignif::geom_signif(
      comparisons = pc$comparisons, annotations = pc$annotations,
      map_signif_level = FALSE, tip_length = 0.015,
      textsize = 4, color = "grey30", step_increase = 0.10)
  p
}

fig4b_mox <- make_phage_box(
  filter(phage_tpm_df, PhageType == "Methyloprofundus phage"),
  dt_mox, kw_mox$p.value,
  "Methyloprofundus Phage", "#C62828")
fig4b_cam <- make_phage_box(
  filter(phage_tpm_df, PhageType == "Campylobacter.(novel) phage"),
  dt_cam, kw_cam$p.value,
  "Campylobacter.(novel) Phage", "#0072B2")
fig4b <- patchwork::wrap_plots(fig4b_mox, fig4b_cam, nrow = 1)

# ── Panel C: pooled MOX coupling scatter ────────────────────
mox_mag_id <- mag_taxonomy %>%
  filter(genus_label == "Methyloprofundus") %>% pull(MAG_ID)

mox_mag_per_sample <- tibble(
  SampleID = colnames(mag_tpm_mat),
  MOX_MAG_TPM = as.numeric(mag_tpm_mat[mox_mag_id, ])) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID")

mox_phage_per_sample <- tibble(
  SampleID = colnames(tpm_mat),
  MOX_phage_TPM = colSums(tpm_mat[mox_ids, , drop = FALSE])) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID")

coupling_df <- mox_phage_per_sample %>%
  left_join(mox_mag_per_sample %>% dplyr::select(SampleID, MOX_MAG_TPM),
            by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels))

sp_all <- cor.test(coupling_df$MOX_MAG_TPM, coupling_df$MOX_phage_TPM,
                   method = "spearman")

fig4c <- ggplot(coupling_df,
                aes(log1p(MOX_MAG_TPM), log1p(MOX_phage_TPM), color = Stage)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1),
              color = "grey45", linewidth = 0.8, alpha = 0.15) +
  geom_point(size = 3.2, alpha = 0.85) +
  scale_color_manual(values = stage_pal, name = "Stage") +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Spearman \u03c1 = %.2f\np = %.3f\nn = %d",
                           sp_all$estimate, sp_all$p.value, nrow(coupling_df)),
           hjust = 1.05, vjust = -0.1, size = 3.4,
           fill = "white", color = "grey20", label.size = 0.3) +
  labs(x = "log(1 + Methyloprofundus MAG TPM)",
       y = "log(1 + Methyloprofundus-phage TPM)",
       title = "Virus\u2013Host Density Coupling",
       subtitle = "Methyloprofundus MAG vs cognate phage") +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        legend.position = c(0.12, 0.80),
        legend.background = element_rect(fill = "white", color = "grey80",
                                         linewidth = 0.3))

# ── Panel D: MAG abundance heatmap ──────────────────────────
hm_mag <- log2(mag_tpm_mat[, sample_order_vec] + 0.5)
hm_mag <- t(scale(t(hm_mag)))
hm_mag[hm_mag >  3] <-  3; hm_mag[hm_mag < -3] <- -3

row_ann_mag <- mag_taxonomy %>%
  dplyr::select(MAG_ID, genus_label, phylum) %>%
  column_to_rownames("MAG_ID") %>%
  dplyr::rename(Genus = genus_label, Phylum = phylum)
row_ann_mag <- row_ann_mag[rownames(hm_mag), , drop = FALSE]

phylum_uniq <- unique(row_ann_mag$Phylum)
phylum_cols <- setNames(brewer.pal(max(3, length(phylum_uniq)), "Set1")[seq_along(phylum_uniq)],
                        phylum_uniq)
genus_uniq  <- unique(row_ann_mag$Genus)
others <- setdiff(genus_uniq, c("Methyloprofundus","Campylobacter.(novel)"))
genus_cols  <- setNames(
  c("#C62828", "#0072B2",
    colorRampPalette(brewer.pal(8, "Dark2"))(length(others))),
  c("Methyloprofundus","Campylobacter.(novel)", others))

mag_ann_colors <- list(Stage = stage_pal,
                       Phylum = phylum_cols, Genus = genus_cols)

pdf(file.path(FIG_DIR, "Fig4D_MAG_heatmap.pdf"), width = 12, height = 8)
pheatmap(hm_mag,
         color  = colorRampPalette(c("#2166AC","white","#B2182B"))(101),
         breaks = seq(-3, 3, length.out = 102),
         cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = col_ann, annotation_row = row_ann_mag,
         annotation_colors = mag_ann_colors,
         labels_row = row_ann_mag$Genus,
         show_colnames = FALSE,
         fontsize_row = 9, cellwidth = 6, cellheight = 14,
         border_color = NA,
         gaps_col = cumsum(table(col_ann$Stage)[stage_levels])[1:2],
         main = "Bacterial MAG Abundance (z-score)")
dev.off()

# Combine A-C
fig4_abc <- patchwork::wrap_plots(fig4a, fig4b, fig4c, ncol = 3,
                                  widths = c(1, 1.2, 1.2)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 4  |  Host Prediction & Virus\u2013Host Dynamics",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(FIG_DIR, "Fig4_ABC.pdf"), fig4_abc,
       width = 16, height = 6, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig4_ABC.png"), fig4_abc,
       width = 16, height = 6, dpi = 300)
cat("Figure 4 saved\n")
