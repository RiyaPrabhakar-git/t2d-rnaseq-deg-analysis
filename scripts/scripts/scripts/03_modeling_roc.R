library(pROC)

vsd <- vst(dds)

res <- results(dds)
res_df <- as.data.frame(res)

biomarkers <- res_df %>% arrange(padj) %>% head(5)

# ROC model
