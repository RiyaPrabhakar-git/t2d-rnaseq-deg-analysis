# ==============================
# 04_visualization.R
# PCA, volcano, and heatmap
# ==============================

source("scripts/03_deseq2_analysis.R")

vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, "results/r_objects/vst_object.rds")

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "disease_group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = disease_group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot: T2D vs Non-Diabetic Samples")

ggsave("results/figures/pca_plot.png", pca_plot, width = 8, height = 6)

# Volcano plot
volcano_data <- res_df %>%
  mutate(neg_log10_padj = -log10(padj))

volcano_plot <- ggplot(
  volcano_data,
  aes(x = log2FoldChange, y = neg_log10_padj, color = significance)
) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: T2D vs Non-Diabetic",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

ggsave("results/figures/volcano_plot.png", volcano_plot, width = 8, height = 6)

# Top DEG heatmap
top_genes <- sig_genes %>%
  arrange(padj) %>%
  slice_head(n = 30) %>%
  pull(gene)

heat_mat <- assay(vsd)[top_genes, ]
heat_mat <- heat_mat - rowMeans(heat_mat)

annotation_col <- data.frame(
  disease_group = meta_matched$disease_group
)

rownames(annotation_col) <- meta_matched$sample_id

pheatmap(
  heat_mat,
  annotation_col = annotation_col,
  scale = "row",
  filename = "results/figures/top30_deg_heatmap.png"
)
