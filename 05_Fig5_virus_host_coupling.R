# ============================================================
# 05_Fig5_virus_host_coupling.R
# Figure 5: stage-stratified coupling, multi-host rho heatmap,
#           VMR and LME model
# ============================================================

# ── Multi-host coupling table ───────────────────────────────
mag_by_genus <- mag_tpm_mat %>%
  as.data.frame() %>%
  rownames_to_column("MAG_ID") %>%
  left_join(dplyr::select(mag_taxonomy, MAG_ID, genus_label), by = "MAG_ID") %>%
  pivot_longer(-c(MAG_ID, genus_label),
               names_to = "SampleID", values_to = "MAG_TPM") %>%
  group_by(genus_label, SampleID) %>%
  summarise(MAG_TPM = sum(MAG_TPM), .groups = "drop")

phage_by_genus <- host_lookup %>%
  filter(Contig_norm %in% rownames(tpm_mat)) %>%
  group_by(host_genus) %>%
  group_modify(~ {
    ids <- .x$Contig_norm
    tibble(SampleID = colnames(tpm_mat),
           phage_TPM = colSums(tpm_mat[ids, , drop = FALSE]))
  }) %>%
  ungroup() %>%
  dplyr::rename(genus_label = host_genus)

coupling_multi <- mag_by_genus %>%
  inner_join(phage_by_genus, by = c("genus_label", "SampleID")) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels))

# ── Panel A: MOX stratified scatter ─────────────────────────
mox_df <- coupling_multi %>% filter(genus_label == "Methyloprofundus")

stage_rho_mox <- mox_df %>%
  group_by(Stage) %>%
  summarise(rho = cor(MAG_TPM, phage_TPM, method = "spearman"),
            p   = cor.test(MAG_TPM, phage_TPM, method = "spearman")$p.value,
            n   = n(), .groups = "drop") %>%
  mutate(label = sprintf("\u03c1 = %.2f\np = %.3f", rho, p))

fig5a <- ggplot(mox_df, aes(log1p(MAG_TPM), log1p(phage_TPM))) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#C62828", fill = "#C62828",
              alpha = 0.15, linewidth = 0.8) +
  geom_point(aes(color = Stage), size = 2.8, alpha = 0.85) +
  geom_text(data = stage_rho_mox, aes(label = label),
            x = Inf, y = -Inf, hjust = 1.08, vjust = -0.2,
            size = 3.2, color = "grey25", lineheight = 0.85,
            inherit.aes = FALSE) +
  facet_wrap(~ Stage, nrow = 1) +
  scale_color_manual(values = stage_pal, guide = "none") +
  labs(x = "log(1 + Methyloprofundus MAG TPM)",
       y = "log(1 + Methyloprofundus-phage TPM)",
       title = "Virus\u2013Host Coupling: Methyloprofundus",
       subtitle = "Stage-stratified Spearman correlation") +
  theme_bw(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 11),
        strip.text       = element_text(face = "bold", size = 11, color = "white"),
        strip.background = element_rect(fill = "#C62828"))

# ── Panel B: multi-host rho heatmap ─────────────────────────
rho_mat_df <- coupling_multi %>%
  group_by(genus_label, Stage) %>%
  summarise(
    rho = ifelse(n() >= 5 & sd(MAG_TPM) > 0 & sd(phage_TPM) > 0,
                 cor(MAG_TPM, phage_TPM, method = "spearman"),
                 NA_real_),
    p   = ifelse(n() >= 5 & sd(MAG_TPM) > 0 & sd(phage_TPM) > 0,
                 cor.test(MAG_TPM, phage_TPM, method = "spearman")$p.value,
                 NA_real_),
    .groups = "drop"
  ) %>%
  mutate(sig_lab = sig_label(p))

rho_mat <- rho_mat_df %>%
  dplyr::select(genus_label, Stage, rho) %>%
  pivot_wider(names_from = Stage, values_from = rho) %>%
  column_to_rownames("genus_label") %>% as.matrix()
rho_mat <- rho_mat[, stage_levels, drop = FALSE]

rho_long <- rho_mat_df %>%
  mutate(Stage = factor(Stage, levels = stage_levels),
         genus_label = factor(genus_label,
                              levels = rownames(rho_mat)[
                                order(rowMeans(rho_mat, na.rm = TRUE))]))

fig5b <- ggplot(rho_long, aes(Stage, genus_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_lab), size = 5, vjust = 0.75, color = "grey20") +
  geom_text(aes(label = ifelse(!is.na(rho), sprintf("%.2f", rho), "NA")),
            size = 3.0, vjust = -0.5, color = "grey20") +
  scale_fill_gradientn(
    colors  = c("#2166AC","#92C5DE","white","#F4A582","#B2182B"),
    values  = scales::rescale(c(-1, -0.3, 0, 0.3, 1)),
    limits  = c(-1, 1), na.value = "grey88", name = "Spearman \u03c1") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title    = "Multi-Host Virus\u2013Host Coupling",
       subtitle = "Spearman \u03c1 per stage  (* p<0.05, ** p<0.01, *** p<0.001)") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        axis.text.x   = element_text(face = "bold", size = 10.5),
        panel.grid    = element_blank(),
        legend.key.height = unit(1.2, "cm"))

# ── Panel C: VMR ────────────────────────────────────────────
vmr_df <- tibble(
  SampleID  = colnames(tpm_mat),
  phage_TPM = colSums(tpm_mat),
  mag_TPM   = colSums(mag_tpm_mat[, colnames(tpm_mat)]),
  log_VMR   = log10((phage_TPM + 1) / (mag_TPM + 1))
) %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage), by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels))

kw_vmr <- kruskal.test(log_VMR ~ Stage, data = vmr_df)
dt_vmr <- dunn.test::dunn.test(vmr_df$log_VMR, vmr_df$Stage,
                               method = "BH", list = FALSE, kw = FALSE,
                               label = FALSE)

pc_vmr <- list(
  comparisons = lapply(dt_vmr$comparisons[dt_vmr$P.adjusted < 0.05],
                       function(s) trimws(strsplit(s, " - ")[[1]])),
  annotations = sig_label(dt_vmr$P.adjusted[dt_vmr$P.adjusted < 0.05]))

fig5c <- ggplot(vmr_df, aes(Stage, log_VMR, fill = Stage, color = Stage)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = 0.55, width = 0.52) +
  geom_jitter(width = 0.18, size = 2.0, alpha = 0.75,
              shape = 21, stroke = 0.25, color = "white") +
  {if (length(pc_vmr$comparisons) > 0)
    ggsignif::geom_signif(comparisons = pc_vmr$comparisons,
                          annotations = pc_vmr$annotations,
                          map_signif_level = FALSE, tip_length = 0.015,
                          textsize = 4, color = "grey30", step_increase = 0.10)} +
  annotate("text", x = 2, y = Inf,
           label = sprintf("KW p = %.3f", kw_vmr$p.value),
           vjust = -0.4, size = 3.4, fontface = "italic",
           color = ifelse(kw_vmr$p.value < 0.05, "grey20", "grey55")) +
  scale_fill_manual(values = stage_pal, guide = "none") +
  scale_color_manual(values = stage_pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "log\u2081\u2080(Phage TPM / Bacterial MAG TPM)",
       title = "Virus-to-Microbe Ratio",
       subtitle = "Total phage TPM / total MAG TPM per sample") +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        panel.grid.major.x = element_blank())

# ── Panel D: LME ────────────────────────────────────────────
coupling_lme <- tibble(
  SampleID = colnames(tpm_mat),
  MOX_phage_TPM = colSums(tpm_mat[mox_ids, , drop = FALSE])) %>%
  left_join(tibble(SampleID = colnames(mag_tpm_mat),
                   MOX_MAG_TPM = as.numeric(mag_tpm_mat[mox_mag_id, ])),
            by = "SampleID") %>%
  left_join(meta_sub %>% dplyr::select(SampleID, Stage, PrefixGroup),
            by = "SampleID") %>%
  mutate(Stage = factor(Stage, levels = stage_levels),
         log_phage = log1p(MOX_phage_TPM),
         log_mag   = log1p(MOX_MAG_TPM))

lme_mod <- lme4::lmer(
  log_phage ~ log_mag + Stage + (1 | PrefixGroup),
  data = coupling_lme, REML = FALSE)

lme_tidy <- broom.mixed::tidy(lme_mod, conf.int = TRUE, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = case_when(
      term == "log_mag"     ~ "log(MOX-MAG TPM+1)",
      term == "StageMiddle" ~ "Stage: Middle vs Early",
      term == "StageLate"   ~ "Stage: Late vs Early",
      TRUE ~ term),
    sig = sig_label(p.value),
    term_label = factor(term_label,
                        levels = rev(c("log(MOX-MAG TPM+1)",
                                       "Stage: Middle vs Early",
                                       "Stage: Late vs Early"))))

write_csv(lme_tidy, file.path(TABLE_DIR, "TableS7_LME_coefficients.csv"))

fig5d <- ggplot(lme_tidy,
                aes(estimate, term_label,
                    xmin = conf.low, xmax = conf.high,
                    color = p.value < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey55", linewidth = 0.6) +
  geom_errorbarh(height = 0.2, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("\u03b2=%.2f %s", estimate, sig)),
            nudge_y = 0.28, size = 3.4, color = "grey20") +
  scale_color_manual(values = c("TRUE" = "#C62828", "FALSE" = "grey55"),
                     guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(x = "Coefficient estimate (95% CI)", y = NULL,
       title    = "LME Model: MOX-Phage Drivers",
       subtitle = "log(phage) ~ log(MAG) + Stage + (1|PrefixGroup)") +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        panel.grid.major.y = element_blank())

# ── Combine ─────────────────────────────────────────────────
fig5 <- patchwork::wrap_plots(
  fig5a,
  patchwork::wrap_plots(fig5b, fig5c, fig5d, nrow = 1,
                        widths = c(1, 0.85, 0.85)),
  ncol = 1, heights = c(1, 1.1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 5  |  Virus\u2013Host Density Coupling Across Developmental Stages",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(FIG_DIR, "Fig5_coupling.pdf"), fig5,
       width = 17, height = 11, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig5_coupling.png"), fig5,
       width = 17, height = 11, dpi = 300)
cat("Figure 5 saved\n")
