# ==============================
# 06_pathway_enrichment.R
# GO Biological Process enrichment
# ==============================

source("scripts/05_biomarker_modeling_roc.R")

sig_ensembl <- sig_genes$gene
sig_ensembl_clean <- gsub("\\..*", "", sig_ensembl)

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = sig_ensembl_clean,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

entrez_ids <- na.omit(entrez_ids)

ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

go_results <- as.data.frame(ego)

write.csv(go_results, "results/tables/go_biological_process_enrichment.csv", row.names = FALSE)

png("results/figures/go_enrichment_dotplot.png", width = 1200, height = 900, res = 150)
dotplot(ego, showCategory = 15)
dev.off()

saveRDS(ego, "results/r_objects/go_enrichment_object.rds")

print(head(go_results))
