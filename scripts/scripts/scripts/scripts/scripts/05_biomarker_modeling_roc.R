# ==============================
# 05_biomarker_modeling_roc.R
# Biomarker selection, logistic regression, ROC-AUC
# ==============================

source("scripts/04_visualization.R")

biomarkers <- sig_genes %>%
  arrange(padj) %>%
  slice_head(n = 10)

write.csv(biomarkers, "results/tables/top_biomarkers.csv", row.names = FALSE)

biomarker_genes <- biomarkers$gene[1:3]

model_mat <- as.data.frame(t(assay(vsd)[biomarker_genes, ]))
model_mat$sample <- rownames(model_mat)

model_df <- model_mat %>%
  left_join(meta_matched, by = c("sample" = "sample_id")) %>%
  mutate(outcome = ifelse(disease_group == "T2D", 1, 0))

model_input <- model_df[, c(biomarker_genes, "outcome")]
model_input$outcome <- as.numeric(model_input$outcome)

formula_str <- paste("outcome ~", paste(biomarker_genes, collapse = " + "))

logit_model <- glm(
  as.formula(formula_str),
  data = model_input,
  family = binomial()
)

model_input$predicted_probability <- predict(logit_model, type = "response")

roc_obj <- roc(model_input$outcome, model_input$predicted_probability)
auc_value <- auc(roc_obj)

write.csv(model_input, "results/models/roc_predictions.csv", row.names = FALSE)
write.csv(
  data.frame(AUC = as.numeric(auc_value)),
  "results/models/roc_auc_value.csv",
  row.names = FALSE
)

png("results/figures/roc_auc_plot.png", width = 1000, height = 800, res = 150)
plot(
  roc_obj,
  main = paste("ROC Curve: Top 3 Biomarker Model, AUC =", round(auc_value, 3))
)
dev.off()

# Biomarker-only PCA
biomarker_mat <- assay(vsd)[biomarker_genes, ]
pca <- prcomp(t(biomarker_mat))

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  disease_group = meta_matched$disease_group
)

biomarker_pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = disease_group)) +
  geom_point(size = 4) +
  theme_minimal() +
  ggtitle("Biomarker-Based PCA Separation: T2D vs ND")

ggsave(
  "results/figures/biomarker_pca_plot.png",
  biomarker_pca_plot,
  width = 8,
  height = 6
)

saveRDS(logit_model, "results/r_objects/logistic_model.rds")

print(auc_value)
print(summary(logit_model))
