# ============================================================
# 01_Fig1_vOTU_overview.R
# Figure 1: vOTU catalogue characteristics
# ============================================================

n_votu <- nrow(votu_annot)

# ── Panel A: geNomad virus score ────────────────────────────
score_df <- votu_annot %>% filter(!is.na(virus_score))
n_hi   <- sum(score_df$virus_score >= 0.95)
pct_hi <- n_hi / nrow(score_df) * 100

fig1a <- ggplot(score_df, aes(virus_score)) +
  geom_histogram(bins = 20, fill = "#1565C0", color = "white",
                 alpha = 0.88, linewidth = 0.25) +
  geom_vline(xintercept = 0.95, linetype = "dashed",
             color = "#C62828", linewidth = 0.9) +
  annotate("text", x = 0.945, y = Inf,
           label = sprintf("%.0f%% \u2265 0.95\n(n = %d)", pct_hi, n_hi),
           hjust = 1.05, vjust = 2, color = "#C62828",
           size = 3.6, fontface = "bold") +
  scale_x_continuous(limits = c(0.7, 1.0),
                     breaks = c(0.7, 0.8, 0.9, 0.95, 1.0)) +
  labs(x = "geNomad Virus Score", y = "vOTU count",
       title    = "Viral Identification Confidence",
       subtitle = sprintf("geNomad (n = %d)", nrow(score_df))) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9))

# ── Panel B: genome length × lifecycle ──────────────────────
lc_n <- votu_annot %>%
  dplyr::count(lifecycle) %>%
  mutate(leg_lab = sprintf("%s (n=%d)", lifecycle, n))
lc_pal  <- c(lytic = "#1565C0", lysogenic = "#C62828", Unknown = "#B0BEC5")
lc_labs <- setNames(lc_n$leg_lab, lc_n$lifecycle)

fig1b <- ggplot(votu_annot, aes(contig_length / 1000, fill = lifecycle)) +
  geom_histogram(bins = 35, color = "white", alpha = 0.88, linewidth = 0.2) +
  scale_fill_manual(values = lc_pal, name = "Lifecycle", labels = lc_labs) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 50)) +
  labs(x = "Genome length (kb, log scale)", y = "Count",
       title    = "Genome Length Distribution",
       subtitle = sprintf("n = %d vOTUs", n_votu)) +
  theme_bw(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 11),
        plot.subtitle    = element_text(color = "grey45", size = 9),
        legend.position  = c(0.79, 0.72),
        legend.background = element_rect(fill = "white", color = "grey80",
                                         linewidth = 0.3))

# ── Panel C: class-level taxonomy donut ─────────────────────
class_pal <- c(Caudoviricetes = "#1565C0",
               Megaviricetes  = "#2E7D32",
               Unclassified   = "#CFD8DC")
class_cnt <- votu_annot %>%
  dplyr::count(class, sort = TRUE) %>%
  mutate(class_lab = case_when(
    class == "Caudoviricetes" ~ "Caudoviricetes",
    class == "Megaviricetes"  ~ "Megaviricetes",
    TRUE                      ~ "Unclassified"
  )) %>%
  group_by(class_lab) %>% summarise(n = sum(n), .groups = "drop") %>%
  mutate(pct   = n / sum(n) * 100,
         label = sprintf("%s\n%.0f%%\n(n=%d)", class_lab, pct, n),
         class_lab = factor(class_lab,
                            levels = c("Caudoviricetes","Megaviricetes","Unclassified")))

fig1c <- ggplot(class_cnt, aes(x = 2, y = n, fill = class_lab)) +
  geom_col(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3.3, fontface = "bold", color = "white",
            lineheight = 0.9) +
  coord_polar("y") + xlim(0.3, 2.5) +
  scale_fill_manual(values = class_pal, name = "Viral Class") +
  labs(title = "Taxonomic Composition", subtitle = "Class level (geNomad)") +
  theme_void() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
        plot.subtitle = element_text(hjust = 0.5, color = "grey45", size = 9))

# ── Panel D: predicted bacterial hosts ──────────────────────
host_cnt <- host_lookup %>%
  dplyr::count(host_genus, sort = TRUE) %>%
  mutate(pct = n / sum(n) * 100,
         genus_lab = ifelse(row_number() <= 9, host_genus, "Other")) %>%
  group_by(genus_lab) %>%
  summarise(n = sum(n), pct = sum(pct), .groups = "drop") %>%
  arrange(desc(n)) %>%
  mutate(genus_lab = factor(genus_lab,
                            levels = rev(c(setdiff(genus_lab[order(-n)], "Other"),
                                           "Other"))),
         is_mox = genus_lab == "Methyloprofundus")

fig1d <- ggplot(host_cnt, aes(genus_lab, pct, fill = is_mox)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            hjust = -0.12, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#C62828", "FALSE" = "#78909C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(x = NULL, y = "Host prediction proportion (%)",
       title    = "Predicted Bacterial Hosts",
       subtitle = sprintf("iPHoP \u226590%% confidence (n = %d pairs)",
                          nrow(host_pred))) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        panel.grid.major.y = element_blank())

# ── Panel E: AMG functional modules ─────────────────────────
amg_mod_cnt <- amg_hq %>%
  filter(module != "Other") %>%
  dplyr::count(module, source, sort = FALSE) %>%
  group_by(module) %>% mutate(total = sum(n)) %>% ungroup() %>%
  arrange(total) %>%
  mutate(module = factor(module, levels = unique(module)))

fig1e <- ggplot(amg_mod_cnt,
                aes(module, n, fill = module, alpha = source)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(
    data = amg_mod_cnt %>% group_by(module) %>%
      summarise(total = sum(n), .groups = "drop"),
    aes(module, total, label = total, fill = NULL, alpha = NULL),
    hjust = -0.2, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = mod_pal, guide = "none") +
  scale_alpha_manual(values = c("DRAM-v" = 1, "VIBRANT" = 0.55),
                     name = "Source",
                     labels = c("DRAM-v (rank B)", "VIBRANT (supp.)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_flip() +
  labs(x = NULL, y = "AMG count",
       title    = "Auxiliary Metabolic Genes",
       subtitle = sprintf("DRAM-v rank-B + VIBRANT (n = %d AMGs)",
                          nrow(amg_hq %>% filter(module != "Other")))) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey45", size = 9),
        panel.grid.major.y = element_blank(),
        legend.position   = c(0.75, 0.15),
        legend.background = element_rect(fill = "white", color = "grey80",
                                         linewidth = 0.3))

# ── Combine ─────────────────────────────────────────────────
fig1_top <- patchwork::wrap_plots(fig1a, fig1b, fig1c, nrow = 1)
fig1_bot <- patchwork::wrap_plots(fig1d, fig1e, nrow = 1, widths = c(1.1, 1))

fig1 <- patchwork::wrap_plots(fig1_top, fig1_bot, ncol = 1,
                              heights = c(1, 1.15)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 1  |  Deep-sea Virome Catalog Characteristics",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0,
                                             margin = margin(b = 6))))

ggsave(file.path(FIG_DIR, "Fig1_vOTU_overview.pdf"), fig1,
       width = 16, height = 11, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig1_vOTU_overview.png"), fig1,
       width = 16, height = 11, dpi = 300)
cat("Figure 1 saved\n")
