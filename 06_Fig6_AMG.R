# ============================================================
# 06_Fig6_AMG.R
# Figure 6: AMG functional analysis
# ============================================================

# ── Panel A: AMG modules ────────────────────────────────────
mod_cnt <- amg_hq %>%
  filter(module != "Other") %>%
  dplyr::count(module, source) %>%
  group_by(module) %>% mutate(total = sum(n)) %>% ungroup() %>%
  arrange(total) %>%
  mutate(module = factor(module, levels = unique(module)))

fig6a <- ggplot(mod_cnt, aes(module, n, fill = module, alpha = source)) +
  geom_col(width = 0.70, color = "white", linewidth = 0.3) +
  geom_text(
    data = mod_cnt %>% group_by(module) %>%
      summarise(total = sum(n), .groups = "drop"),
    aes(module, total, label = total, fill = NULL, alpha = NULL),
    hjust = -0.2, size = 3.6, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = mod_pal, guide = "none") +
  scale_alpha_manual(values = c("DRAM-v" = 1, "VIBRANT" = 0.5),
                     name = "Source",
                     labels = c("DRAM-v (rank B)", "VIBRANT (supp.)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_flip() +
  labs(x = NULL, y = "AMG count", title = "AMG Functional Modules") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_blank(),
        legend.position = c(0.72, 0.15),
        legend.background = element_rect(fill = "white", color = "grey80",
                                         linewidth = 0.3))

# ── Panel B: host_symbiont x module bubble ──────────────────
top_sym <- amg %>%
  filter(!is.na(host_symbiont), module != "Other") %>%
  dplyr::count(host_symbiont, sort = TRUE) %>%
  slice_head(n = 8) %>% pull(host_symbiont)

bubble_df <- amg %>%
  filter(module != "Other") %>%
  mutate(sym_lab = case_when(
    host_symbiont %in% top_sym ~ host_symbiont,
    is.na(host_symbiont)       ~ "Unknown",
    TRUE                       ~ "Other bacteria"),
    sym_lab = factor(sym_lab, levels = c(top_sym, "Other bacteria", "Unknown"))) %>%
  dplyr::count(module, sym_lab) %>%
  mutate(module = factor(module, levels = names(mod_pal)))

bath_idx <- which(levels(bubble_df$sym_lab) == "Bathymodiolus symbiont")

fig6b <- ggplot(bubble_df, aes(sym_lab, module, size = n, fill = module)) +
  geom_point(shape = 21, alpha = 0.85, color = "white", stroke = 0.3) +
  geom_text(aes(label = ifelse(n >= 5, n, "")),
            size = 2.8, color = "white", fontface = "bold") +
  {if (length(bath_idx) > 0)
    annotate("rect", xmin = bath_idx - 0.5, xmax = bath_idx + 0.5,
             ymin = 0.4, ymax = length(unique(bubble_df$module)) + 0.6,
             fill = NA, color = "#C62828", linewidth = 1.0)} +
  scale_fill_manual(values = mod_pal, guide = "none") +
  scale_size_continuous(name = "AMG count", range = c(2, 12),
                        breaks = c(5, 20, 50, 100)) +
  labs(x = NULL, y = NULL,
       title    = "AMG Distribution by Predicted Host Symbiont",
       subtitle = "Red box = Bathymodiolus symbiont") +
  theme_bw(base_size = 10.5) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "#C62828", size = 9, face = "italic"),
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 9))

# ── Panel C: Bathymodiolus-symbiont AMG detail ──────────────
bath_amg <- amg %>%
  filter(host_symbiont == "Bathymodiolus symbiont", module != "Other") %>%
  dplyr::count(module, ko_id, ko_name) %>%
  mutate(ko_short = str_extract(ko_name, "^[^;\\[]+") %>%
         str_trim() %>% str_sub(1, 35),
         label = sprintf("%s\n(%s)", ko_short, ko_id),
         module = factor(module, levels = names(mod_pal))) %>%
  arrange(module, desc(n)) %>%
  mutate(label = factor(label, levels = rev(unique(label))))

fig6c <- ggplot(bath_amg, aes(label, n, fill = module)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.4, fontface = "bold") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = mod_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  coord_flip() +
  labs(x = NULL, y = "AMG copy number",
       title    = "Bathymodiolus Symbiont-Targeting AMGs",
       subtitle = sprintf("n = %d AMGs", sum(bath_amg$n))) +
  theme_bw(base_size = 10) +
  theme(plot.title       = element_text(face = "bold", size = 11),
        strip.text       = element_text(face = "bold", size = 8.5, color = "white"),
        strip.background = element_rect(fill = "#37474F"),
        panel.grid.major.y = element_blank())

# ── Combine A-C; D heatmap saved separately ─────────────────
fig6_abc <- patchwork::wrap_plots(fig6a, fig6b, fig6c, ncol = 3,
                                  widths = c(0.85, 1.3, 0.85)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 6  |  Auxiliary Metabolic Genes in the Deep-Sea Mussel Phageome",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(FIG_DIR, "Fig6_ABC.pdf"), fig6_abc,
       width = 18, height = 7, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig6_ABC.png"), fig6_abc,
       width = 18, height = 7, dpi = 300)

# ── Panel D: AMG-carrying vOTU heatmap ──────────────────────
amg_votus <- amg_hq %>%
  filter(module != "Other", scaffold %in% rownames(tpm_mat)) %>%
  dplyr::count(scaffold, sort = TRUE) %>%
  slice_head(n = 40) %>% pull(scaffold)

if (length(amg_votus) >= 5) {
  hm_amg <- tpm_mat[amg_votus, sample_order_vec, drop = FALSE]
  hm_amg <- log2(hm_amg + 0.5); hm_amg <- t(scale(t(hm_amg)))
  hm_amg[hm_amg >  3] <-  3; hm_amg[hm_amg < -3] <- -3

  row_ann_amg <- amg_hq %>%
    filter(scaffold %in% amg_votus) %>%
    group_by(scaffold) %>%
    summarise(Top_module = names(sort(table(module), decreasing = TRUE))[1],
              Bath_symbiont = any(host_symbiont == "Bathymodiolus symbiont", na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Host_type = if_else(Bath_symbiont,
                               "Bathymodiolus symbiont", "Other host")) %>%
    column_to_rownames("scaffold") %>%
    dplyr::select(Top_module, Host_type)

  amg_ko_lab <- amg_hq %>%
    filter(scaffold %in% rownames(row_ann_amg)) %>%
    group_by(scaffold) %>%
    summarise(kos = paste(unique(ko_id)[1:min(2, n_distinct(ko_id))], collapse = "/"),
              .groups = "drop") %>% deframe()
  row_labels_amg <- sprintf("%s | %s",
                            amg_ko_lab[rownames(row_ann_amg)],
                            row_ann_amg$Host_type %>%
                              str_replace("Bathymodiolus symbiont", "Bath.sym") %>%
                              str_replace("Other host", "Other"))
  names(row_labels_amg) <- rownames(row_ann_amg)

  mod_present  <- intersect(unique(row_ann_amg$Top_module), names(mod_pal))
  missing_mods <- setdiff(unique(row_ann_amg$Top_module), names(mod_pal))
  mod_cols_amg <- mod_pal[mod_present]
  if (length(missing_mods) > 0)
    mod_cols_amg <- c(mod_cols_amg,
                      setNames(rep("#90A4AE", length(missing_mods)), missing_mods))

  amg_ann_colors <- list(
    Stage = stage_pal,
    Top_module = mod_cols_amg,
    Host_type = c("Bathymodiolus symbiont" = "#C62828",
                  "Other host" = "#78909C"))

  pdf(file.path(FIG_DIR, "Fig6D_AMG_heatmap.pdf"), width = 12, height = 10)
  pheatmap(hm_amg[rownames(row_ann_amg), , drop = FALSE],
           color  = colorRampPalette(c("#2166AC","white","#B2182B"))(101),
           breaks = seq(-3, 3, length.out = 102),
           cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = col_ann, annotation_row = row_ann_amg,
           annotation_colors = amg_ann_colors,
           labels_row = row_labels_amg[rownames(row_ann_amg)],
           show_colnames = FALSE,
           fontsize_row = 7.5, cellwidth = 6, cellheight = 10,
           border_color = NA,
           gaps_col = cumsum(table(col_ann$Stage)[stage_levels])[1:2],
           main = "AMG-carrying vOTU Abundance (z-score)")
  dev.off()
}

cat("Figure 6 saved\n")
