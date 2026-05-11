# ============================================================
# 02_Fig2_diversity.R
# Figure 2: viral community diversity across stages
# ============================================================

# ── Community matrix ────────────────────────────────────────
comm <- t(tpm_mat)
set.seed(42)
bc_dist <- vegan::vegdist(comm, method = "bray")

# PERMANOVA & ANOSIM
perm_res   <- vegan::adonis2(bc_dist ~ Stage, data = meta_sub, permutations = 999)
anosim_res <- vegan::anosim(bc_dist, meta_sub$Stage, permutations = 999)
perm_R2 <- perm_res$R2[1]; perm_p <- perm_res$`Pr(>F)`[1]

# ── PCoA ────────────────────────────────────────────────────
pcoa_res <- ape::pcoa(bc_dist)
var_exp  <- pcoa_res$values$Relative_eig[1:2] * 100
pcoa_df  <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df  <- pcoa_df %>%
  rownames_to_column("SampleID") %>%
  left_join(meta_sub, by = "SampleID")

# Outlier handling for axis truncation
pc2_q   <- quantile(pcoa_df$PC2, c(0.25, 0.75))
thresh  <- pc2_q[1] - 3 * (pc2_q[2] - pc2_q[1])
pcoa_main <- pcoa_df %>% filter(PC2 >= thresh)
n_out  <- sum(pcoa_df$PC2 < thresh)

xpad  <- diff(range(pcoa_main$PC1)) * 0.10
ypad  <- diff(range(pcoa_main$PC2)) * 0.12
xlims <- range(pcoa_main$PC1) + c(-xpad, xpad)
ylims <- range(pcoa_main$PC2) + c(-ypad, ypad)

stage_shape <- c(Early = 21, Middle = 24, Late = 22)
stats_lab   <- sprintf(
  "PERMANOVA: R\u00b2 = %.3f, p = %.3f\nANOSIM: R = %.3f, p = %.3f",
  perm_R2, perm_p, anosim_res$statistic, anosim_res$signif)

fig2a <- ggplot(pcoa_df, aes(PC1, PC2, fill = Stage, shape = Stage)) +
  stat_ellipse(aes(color = Stage), type = "t", level = 0.95,
               linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  geom_point(size = 3.2, alpha = 0.88, stroke = 0.3, color = "white") +
  {if (n_out > 0)
    annotate("text", x = xlims[2], y = ylims[1],
             label = sprintf("%d outlier%s outside view",
                             n_out, ifelse(n_out > 1, "s", "")),
             hjust = 1, vjust = 0, size = 3, color = "grey50",
             fontface = "italic")} +
  annotate("label", x = xlims[2], y = ylims[2], label = stats_lab,
           hjust = 1, vjust = 1, size = 3.2, lineheight = 1.2,
           fill = "white", color = "grey20", label.size = 0.35) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  scale_fill_manual(values = stage_pal, name = "Stage") +
  scale_color_manual(values = stage_pal, name = "Stage") +
  scale_shape_manual(values = stage_shape, name = "Stage") +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       title    = "Viral Community Beta-Diversity",
       subtitle = "PCoA of Bray-Curtis dissimilarity") +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        legend.position = c(0.08, 0.20),
        legend.background = element_rect(fill = "white", color = "grey80",
                                         linewidth = 0.3))

# ── Alpha diversity ─────────────────────────────────────────
alpha_df <- tibble(
  SampleID = colnames(tpm_mat),
  Shannon  = vegan::diversity(t(tpm_mat), "shannon"),
  Richness = colSums(tpm_mat > 0),
  Evenness = vegan::diversity(t(tpm_mat), "shannon") / log(colSums(tpm_mat > 0))
) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels))

write_csv(alpha_df, file.path(TABLE_DIR, "TableS4_alpha_diversity.csv"))
write_csv(broom::tidy(perm_res), file.path(TABLE_DIR, "PERMANOVA_result.csv"))

# ── Boxplot helper ──────────────────────────────────────────
parse_comps <- function(dt, alpha = 0.05) {
  keep <- dt$P.adjusted < alpha
  list(
    comparisons = lapply(dt$comparisons[keep],
                         function(s) trimws(strsplit(s, " - ")[[1]])),
    annotations = sig_label(dt$P.adjusted[keep])
  )
}

make_box <- function(df, yvar, ylabel, kw_p, dt_obj) {
  pc <- parse_comps(dt_obj)
  p <- ggplot(df, aes(Stage, .data[[yvar]],
                      fill = Stage, color = Stage)) +
    geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = 0.55, width = 0.52) +
    geom_jitter(width = 0.18, size = 2.0, alpha = 0.75,
                shape = 21, stroke = 0.25, color = "white") +
    scale_fill_manual(values = stage_pal, guide = "none") +
    scale_color_manual(values = stage_pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    annotate("text", x = 2, y = Inf,
             label = sprintf("KW p = %.3f%s", kw_p,
                             ifelse(kw_p < 0.05, "", " (ns)")),
             vjust = -0.3, size = 3.4, fontface = "italic",
             color = ifelse(kw_p < 0.05, "grey20", "grey55")) +
    labs(x = NULL, y = ylabel) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
  if (length(pc$comparisons) > 0)
    p <- p + ggsignif::geom_signif(
      comparisons = pc$comparisons, annotations = pc$annotations,
      map_signif_level = FALSE, tip_length = 0.015,
      textsize = 4, color = "grey30", step_increase = 0.10)
  p
}

kw_shan <- kruskal.test(Shannon  ~ Stage, data = alpha_df)
kw_rich <- kruskal.test(Richness ~ Stage, data = alpha_df)
kw_even <- kruskal.test(Evenness ~ Stage, data = alpha_df)

dt_shan <- dunn.test::dunn.test(alpha_df$Shannon,  alpha_df$Stage,
                                method = "BH", list = FALSE, kw = FALSE,
                                label = FALSE, wrap = FALSE)
dt_rich <- dunn.test::dunn.test(alpha_df$Richness, alpha_df$Stage,
                                method = "BH", list = FALSE, kw = FALSE,
                                label = FALSE, wrap = FALSE)
dt_even <- dunn.test::dunn.test(alpha_df$Evenness, alpha_df$Stage,
                                method = "BH", list = FALSE, kw = FALSE,
                                label = FALSE, wrap = FALSE)

fig2b <- make_box(alpha_df, "Shannon",  "Shannon Index",     kw_shan$p.value, dt_shan)
fig2c <- make_box(alpha_df, "Richness", "Viral Richness",    kw_rich$p.value, dt_rich)
fig2d <- make_box(alpha_df, "Evenness", "Pielou's Evenness", kw_even$p.value, dt_even)

fig2 <- patchwork::wrap_plots(
  fig2a, patchwork::wrap_plots(fig2b, fig2c, fig2d, nrow = 1),
  ncol = 1, heights = c(1.6, 1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 2  |  Viral Community Diversity Across Developmental Stages",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(FIG_DIR, "Fig2_diversity.pdf"), fig2,
       width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig2_diversity.png"), fig2,
       width = 14, height = 10, dpi = 300)
cat("Figure 2 saved\n")
