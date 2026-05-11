# ============================================================
# 09_graphical_abstract.R
# ============================================================

catalog_bar <- tibble(
  category = factor(
    c("Total vOTUs","High/Med quality","With host pred","AMG-carrying"),
    levels = c("Total vOTUs","High/Med quality","With host pred","AMG-carrying")),
  n = c(nrow(tpm_mat),
        sum(votu_annot$checkv_quality %in% c("High-quality","Medium-quality")),
        nrow(host_lookup %>% distinct(Contig_norm)),
        n_distinct(amg$scaffold)))

ga_bar <- ggplot(catalog_bar, aes(category, n, fill = category)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(n)),
            vjust = -0.4, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = c("#1565C0","#0288D1","#00838F","#C62828")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::label_comma()) +
  scale_x_discrete(labels = c("Total\nvOTUs","High/Med\nquality",
                               "With host\nprediction","AMG-\ncarrying")) +
  labs(x = NULL, y = "Count", title = "Deep-Sea Mussel Virome Catalog") +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major.x = element_blank())

ga_donut <- ggplot(host_cnt_full, aes(x = 2, y = n, fill = genus_lab)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct > 5,
                               sprintf("%s\n%.0f%%", genus_lab, pct), "")),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "white", fontface = "bold", lineheight = 0.85) +
  coord_polar("y") + xlim(0.5, 2.5) +
  scale_fill_manual(values = donut_cols) +
  labs(title = "Predicted Bacterial Hosts") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5),
        legend.position = "none")

amg_ga <- amg_hq %>%
  filter(module != "Other") %>%
  dplyr::count(module, sort = TRUE) %>%
  slice_head(n = 6) %>%
  mutate(module = str_wrap(module, 18),
         module = factor(module, levels = rev(unique(module))))

ga_amg <- ggplot(amg_ga, aes(module, n, fill = module)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_flip() +
  labs(x = NULL, y = "AMG count", title = "Top AMG Modules") +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank())

ga_scatter <- ggplot(coupling_lme,
                     aes(log_mag, log_phage, color = Stage)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1),
              color = "#C62828", alpha = 0.15, linewidth = 0.9) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = stage_pal, name = NULL) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("\u03c1 = %.2f", sp_all$estimate),
           hjust = -0.1, vjust = 1.3, size = 4.5, fontface = "bold",
           color = "#C62828") +
  labs(x = "log(Methyloprofundus MAG)", y = "log(Phage TPM)",
       title = "Virus\u2013Host Coupling") +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = c(0.75, 0.25))

findings_text <- paste0(
  "\u2022 ", scales::comma(nrow(tpm_mat)), " vOTUs identified\n",
  "  (3 developmental stages, n=48)\n\n",
  "\u2022 Viral diversity \u2193 Early \u2192 Late\n",
  "  (Shannon & Evenness, p<0.05)\n\n",
  "\u2022 Methyloprofundus phage dominant\n\n",
  "\u2022 Virus\u2013host density coupling\n",
  "  (Spearman \u03c1 = ", round(sp_all$estimate, 2),
  ", p = ", format.pval(sp_all$p.value, digits = 2), ")\n\n",
  "\u2022 AMGs target Bathymodiolus\n",
  "  chemosynthetic symbionts\n",
  "  (n=", sum(amg$host_symbiont == "Bathymodiolus symbiont", na.rm = TRUE),
  " AMGs)")

ga_text <- ggplot() +
  annotate("label", x = 0.5, y = 0.5, label = findings_text,
           hjust = 0.5, vjust = 0.5, size = 3.5, lineheight = 1.3,
           fill = "#E3F2FD", color = "#1565C0",
           label.size = 1.2) +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Key Findings") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5,
                                  size = 10.5, color = "#1565C0"))

ga_top <- patchwork::wrap_plots(ga_bar, ga_donut, ga_amg, nrow = 1,
                                widths = c(1, 0.9, 0.9))
ga_bot <- patchwork::wrap_plots(ga_scatter, ga_text, nrow = 1,
                                widths = c(1.2, 1))

graphical_abstract <- patchwork::wrap_plots(ga_top, ga_bot, ncol = 1,
                                            heights = c(1, 1.1)) +
  patchwork::plot_annotation(
    title = "Viral Ecology of the Deep-Sea Mussel Bathymodiolus",
    subtitle = "Stage-dependent virome dynamics, virus\u2013host coupling and AMGs",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")))

ggsave(file.path(FIG_DIR, "Graphical_Abstract.pdf"), graphical_abstract,
       width = 14, height = 9, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Graphical_Abstract.png"), graphical_abstract,
       width = 14, height = 9, dpi = 300)

cat("Graphical Abstract saved\n")
