# ==============================
# 03_deseq2_analysis.R
# Differential gene expression
# ==============================

source("scripts/02_counts_processing.R")

meta_matched$disease_group <- factor(
  meta_matched$disease_group,
  levels = c("ND", "T2D")
)

rownames(meta_matched) <- meta_matched$sample_id

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix_filt),
  colData = meta_matched,
  design = ~ disease_group
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("disease_group", "T2D", "ND"))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df %>%
  mutate(
    significance = case_when(
      !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
      TRUE ~ "Non-significant"
    ),
    regulation = case_when(
      significance == "Significant" & log2FoldChange > 1 ~ "Upregulated",
      significance == "Significant" & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "No major change"
    )
  )

# Add gene symbols from Ensembl IDs
res_df$ensembl_id <- gsub("\\..*", "", res_df$gene)

res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

write.csv(res_df, "results/tables/deseq2_results_with_symbols.csv", row.names = FALSE)

sig_genes <- res_df %>%
  filter(significance == "Significant")

non_sig_genes <- res_df %>%
  filter(significance == "Non-significant")

write.csv(sig_genes, "results/tables/significant_genes.csv", row.names = FALSE)
write.csv(non_sig_genes, "results/tables/non_significant_genes.csv", row.names = FALSE)

saveRDS(dds, "results/r_objects/dds_object.rds")

print(table(res_df$significance))
print(table(res_df$regulation))
