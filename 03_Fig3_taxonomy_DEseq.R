# ============================================================
# 03_Fig3_taxonomy_DEseq.R
# Figure 3: viral family composition + DESeq2 differential abundance
# ============================================================

# ── TPM long table with family annotation ───────────────────
tpm_long <- tpm_raw %>%
  pivot_longer(-Contig_norm, names_to = "SampleID", values_to = "TPM") %>%
  filter(SampleID %in% meta_sub$SampleID) %>%
  left_join(genomad_tax, by = "Contig_norm") %>%
  mutate(family_label = replace_na(family_label, "Unclassified"))

top8_fam <- tpm_long %>%
  filter(!family_label %in% c("Unclassified","Unclassified Caudoviricetes")) %>%
  group_by(family_label) %>% summarise(tot = sum(TPM), .groups = "drop") %>%
  arrange(desc(tot)) %>% slice_head(n = 8) %>% pull(family_label)

fam_levels <- c(rev(top8_fam), "Unclassified Caudoviricetes", "Unclassified")
fam_colors <- setNames(
  c(brewer.pal(8, "Paired"), "#78909C", "#CFD8DC"),
  fam_levels)

fam_prop <- tpm_long %>%
  mutate(fam_plot = case_when(
    family_label %in% top8_fam                    ~ family_label,
    family_label == "Unclassified Caudoviricetes" ~ "Unclassified Caudoviricetes",
    TRUE                                          ~ "Unclassified")) %>%
  group_by(SampleID, fam_plot) %>%
  summarise(tot = sum(TPM), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(prop = tot / sum(tot) * 100) %>% ungroup() %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels),
         fam_plot = factor(fam_plot, levels = fam_levels))

# ── Panel A: per-sample composition ─────────────────────────
fig3a <- ggplot(
  fam_prop %>% mutate(SampleID = factor(SampleID, levels = sample_order_vec)),
  aes(SampleID, prop, fill = fam_plot)) +
  geom_col(width = 0.92) +
  facet_grid(~ Stage, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fam_colors, name = "Viral Family",
                    guide = guide_legend(ncol = 1, reverse = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative Abundance (% TPM)",
       title = "Sample-Level Viral Family Composition") +
  theme_bw(base_size = 10) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 6.5),
        strip.text       = element_text(face = "bold", size = 11, color = "white"),
        strip.background = element_rect(fill = "#37474F"),
        plot.title       = element_text(face = "bold", size = 11),
        legend.key.size  = unit(0.42, "lines"),
        panel.grid.major.x = element_blank())

# ── Panel B: stage-averaged composition ─────────────────────
fam_stage <- fam_prop %>%
  group_by(Stage, fam_plot) %>%
  summarise(mean_prop = mean(prop),
            se_prop   = sd(prop) / sqrt(n()), .groups = "drop")
total_se <- fam_stage %>%
  group_by(Stage) %>%
  summarise(total_mean = sum(mean_prop),
            total_se   = sqrt(sum(se_prop^2)), .groups = "drop")

fig3b <- ggplot(fam_stage, aes(Stage, mean_prop, fill = fam_plot)) +
  geom_col(position = "stack", width = 0.60, color = "white", linewidth = 0.3) +
  geom_errorbar(data = total_se, inherit.aes = FALSE,
                aes(x = Stage, ymin = total_mean - total_se,
                    ymax = total_mean + total_se),
                width = 0.18, linewidth = 0.7) +
  scale_fill_manual(values = fam_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, 106)) +
  labs(x = NULL, y = "Mean Relative Abundance (%)",
       title = "Stage-Averaged Composition", subtitle = "Error bars = \u00b11 SE") +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        panel.grid.major.x = element_blank())

# ── DESeq2 ──────────────────────────────────────────────────
count_mat <- round(tpm_mat)
count_mat <- count_mat[rowSums(count_mat) > 0, meta_sub$SampleID]
col_data  <- meta_sub %>% column_to_rownames("SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels))

suppressMessages({
  dds  <- DESeq2::DESeqDataSetFromMatrix(count_mat, col_data, ~ Stage)
  keep <- rowSums(DESeq2::counts(dds) >= 5) >= 3
  dds  <- dds[keep, ]
  dds  <- DESeq2::DESeq(dds, quiet = TRUE)
})

mean_tpm <- rowMeans(tpm_mat[rownames(dds), , drop = FALSE])

contrasts_list <- list(
  "Early vs Middle" = c("Stage", "Early", "Middle"),
  "Early vs Late"   = c("Stage", "Early", "Late"),
  "Middle vs Late"  = c("Stage", "Middle", "Late"))

deseq_all <- purrr::imap_dfr(contrasts_list, function(ct, nm) {
  DESeq2::results(dds, contrast = ct, alpha = 0.05, pAdjustMethod = "BH") %>%
    as.data.frame() %>%
    rownames_to_column("Contig_norm") %>%
    filter(!is.na(padj)) %>%
    mutate(comparison = nm,
           mean_tpm = mean_tpm[Contig_norm],
           sig = case_when(
             padj < 0.05 & log2FoldChange >  1 ~ "Up in B",
             padj < 0.05 & log2FoldChange < -1 ~ "Up in A",
             TRUE ~ "NS"))
})

write_csv(deseq_all %>% arrange(comparison, padj),
          file.path(TABLE_DIR, "DESeq2_vOTU_allContrasts.csv"))

# ── Volcano helper ──────────────────────────────────────────
make_volcano <- function(df, title_str, label_a, label_b) {
  n_a <- sum(df$sig == "Up in A"); n_b <- sum(df$sig == "Up in B")
  pal <- c("Up in A" = "#0072B2", "Up in B" = "#D55E00", NS = "grey78")
  ggplot(df, aes(log2FoldChange, -log10(padj),
                 color = sig, size = log10(mean_tpm + 1))) +
    geom_point(alpha = 0.65) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55", linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("\u2191 %s\nn = %d", label_a, n_a),
             hjust = -0.08, vjust = 1.5, color = "#0072B2",
             size = 3.5, fontface = "bold", lineheight = 0.9) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("\u2191 %s\nn = %d", label_b, n_b),
             hjust = 1.08, vjust = 1.5, color = "#D55E00",
             size = 3.5, fontface = "bold", lineheight = 0.9) +
    scale_color_manual(values = pal, guide = "none") +
    scale_size_continuous(name = "log10(mean TPM+1)", range = c(0.5, 3.5),
                          breaks = c(1, 2, 3),
                          guide = guide_legend(override.aes = list(alpha = 0.8))) +
    labs(x = "log\u2082 Fold Change", y = "-log\u2081\u2080(adj. p)", title = title_str) +
    theme_bw(base_size = 10.5) +
    theme(plot.title = element_text(face = "bold", size = 10.5),
          legend.position = "bottom", legend.key.size = unit(0.4, "lines"))
}

fig3d <- make_volcano(filter(deseq_all, comparison == "Early vs Middle"),
                      "Early vs Middle", "Early", "Middle")
fig3e <- make_volcano(filter(deseq_all, comparison == "Early vs Late"),
                      "Early vs Late", "Early", "Late")
fig3f <- make_volcano(filter(deseq_all, comparison == "Middle vs Late"),
                      "Middle vs Late", "Middle", "Late")

# ── Top-40 differential vOTU heatmap (panel C) ──────────────
top40 <- deseq_all %>%
  filter(sig != "NS") %>%
  group_by(Contig_norm) %>%
  summarise(min_padj = min(padj), max_lfc = max(abs(log2FoldChange)),
            .groups = "drop") %>%
  arrange(min_padj, desc(max_lfc)) %>%
  slice_head(n = 40) %>% pull(Contig_norm)

hm_mat <- tpm_mat[top40, sample_order_vec, drop = FALSE]
hm_mat <- log2(hm_mat + 0.5); hm_mat <- t(scale(t(hm_mat)))
hm_mat[hm_mat >  3] <-  3; hm_mat[hm_mat < -3] <- -3

# Top hosts for row annotation
top_hosts <- host_lookup %>%
  filter(host_genus != "Unknown") %>%
  dplyr::count(host_genus, sort = TRUE) %>%
  slice_head(n = 6) %>% pull(host_genus)

host_cols <- setNames(
  c("#C62828",                                # Methyloprofundus
    brewer.pal(6, "Paired")[2:6],              # next 5
    "grey72"),                                  # Unknown
  c(top_hosts, "Unknown"))

row_ann <- data.frame(Contig_norm = top40) %>%
  left_join(host_lookup, by = "Contig_norm") %>%
  left_join(genomad_tax, by = "Contig_norm") %>%
  mutate(
    host_genus   = replace_na(as.character(host_genus),   "Unknown"),
    family_label = replace_na(as.character(family_label), "Unclassified"),
    Host_genus   = if_else(host_genus %in% top_hosts, host_genus, "Unknown"),
    Viral_family = case_when(
      !family_label %in% c("Unclassified","Unclassified Caudoviricetes") ~ family_label,
      TRUE ~ "Unclassified")) %>%
  column_to_rownames("Contig_norm") %>%
  dplyr::select(Host_genus, Viral_family)

named_fams   <- setdiff(unique(row_ann$Viral_family), "Unclassified")
fam_row_cols <- setNames(
  c(brewer.pal(max(3, length(named_fams)), "Set2")[seq_along(named_fams)], "#CFD8DC"),
  c(named_fams, "Unclassified"))

row_labels <- setNames(
  if_else(row_ann$Host_genus != "Unknown",
          sprintf("%s phage", row_ann$Host_genus),
          sprintf("vOTU_%s", str_extract(top40, "\\d+$"))),
  top40)

ann_colors <- list(Stage = stage_pal,
                   Host_genus = host_cols,
                   Viral_family = fam_row_cols)

pdf(file.path(FIG_DIR, "Fig3C_heatmap_top40.pdf"), width = 10, height = 12)
pheatmap(hm_mat,
         color  = colorRampPalette(c("#2166AC","white","#B2182B"))(101),
         breaks = seq(-3, 3, length.out = 102),
         cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = col_ann, annotation_row = row_ann,
         annotation_colors = ann_colors,
         labels_row = row_labels[rownames(hm_mat)],
         show_colnames = FALSE,
         fontsize_row = 7.5, cellwidth = 7, cellheight = 10,
         border_color = NA,
         gaps_col = cumsum(table(col_ann$Stage)[stage_levels])[1:2],
         main = "Top 40 Differential vOTUs (z-score)")
dev.off()

# Combine A/B + volcano panels
fig3_top <- patchwork::wrap_plots(fig3a, fig3b, ncol = 1, heights = c(2.2, 1))
fig3_bot <- patchwork::wrap_plots(fig3d, fig3e, fig3f, nrow = 1)

fig3_panel <- patchwork::wrap_plots(fig3_top, fig3_bot, ncol = 1,
                                    heights = c(1.4, 1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 3  |  Viral Taxonomic Composition and Differential Abundance",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(FIG_DIR, "Fig3_taxonomy.pdf"), fig3_panel,
       width = 16, height = 12, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig3_taxonomy.png"), fig3_panel,
       width = 16, height = 12, dpi = 300)
cat("Figure 3 saved (panel C heatmap saved separately)\n")
