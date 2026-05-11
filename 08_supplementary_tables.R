# ============================================================
# 08_supplementary_tables.R
# Generate supplementary tables S1-S8
# ============================================================

# TableS1 vOTU catalogue
votu_cols <- intersect(
  c("contig_id","contig_length","checkv_quality","completeness","contamination",
    "topology","lifecycle","family_label","virus_score","viral_genes","host_genes",
    "n_hallmarks","provirus"),
  colnames(votu_annot))

tableS1 <- votu_annot %>% dplyr::select(all_of(votu_cols)) %>%
  arrange(desc(factor(checkv_quality,
                      levels = c("High-quality","Medium-quality",
                                 "Low-quality","Not-determined"))),
          desc(contig_length))
write_csv(tableS1, file.path(TABLE_DIR, "TableS1_vOTU_catalogue.csv"))

# TableS2 host prediction
tableS2 <- host_pred %>%
  dplyr::count(host_genus, sort = TRUE) %>%
  mutate(Percentage = round(n / sum(n) * 100, 2),
         Rank = row_number()) %>%
  dplyr::rename(Host_genus = host_genus, Count = n) %>%
  dplyr::select(Rank, Host_genus, Count, Percentage)
write_csv(tableS2, file.path(TABLE_DIR, "TableS2_host_prediction_summary.csv"))

# TableS3 AMG full
tableS3 <- amg %>%
  left_join(host_lookup %>% dplyr::rename(scaffold = Contig_norm),
            by = "scaffold") %>%
  arrange(module, ko_id)
write_csv(tableS3, file.path(TABLE_DIR, "TableS3_AMG_full.csv"))

# TableS4: alpha diversity (already saved by Fig2 script)

# TableS5 beta diversity stats
tableS5 <- bind_rows(
  perm_res %>% as.data.frame() %>%
    rownames_to_column("term") %>% mutate(test = "PERMANOVA"),
  tibble(term = "Stage", statistic = anosim_res$statistic,
         p.value = anosim_res$signif, test = "ANOSIM"))
write_csv(tableS5, file.path(TABLE_DIR, "TableS5_beta_diversity_stats.csv"))

# TableS6: DESeq2 already saved as DESeq2_vOTU_allContrasts.csv → rename
file.copy(file.path(TABLE_DIR, "DESeq2_vOTU_allContrasts.csv"),
          file.path(TABLE_DIR, "TableS6_DESeq2_all_contrasts.csv"),
          overwrite = TRUE)

# TableS7 already saved by Fig5 script

# TableS8 MAG taxonomy + TPM
tableS8 <- mag_taxonomy %>%
  left_join(
    mag_tpm_long %>%
      group_by(MAG_ID, Stage) %>%
      summarise(mean_TPM = mean(MAG_TPM), .groups = "drop") %>%
      pivot_wider(names_from = Stage, values_from = mean_TPM,
                  names_prefix = "mean_TPM_"),
    by = "MAG_ID")
write_csv(tableS8, file.path(TABLE_DIR, "TableS8_MAG_taxonomy_TPM.csv"))

cat("Supplementary tables saved\n")
