# ==============================
# 02_counts_processing.R
# Load counts and align with metadata
# ==============================

source("scripts/01_metadata_processing.R")

counts <- read.delim(
  "data/raw/GSE164416_DP_htseq_counts.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

gene_col <- colnames(counts)[1]

count_matrix <- counts[, -1]
rownames(count_matrix) <- counts[[gene_col]]

count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"

common_samples <- intersect(colnames(count_matrix), meta_clean$sample_id)

count_matrix_sub <- count_matrix[, common_samples]

meta_matched <- meta_clean[match(common_samples, meta_clean$sample_id), ]
meta_matched <- meta_matched[match(colnames(count_matrix_sub), meta_matched$sample_id), ]

stopifnot(all(colnames(count_matrix_sub) == meta_matched$sample_id))

keep <- rowSums(count_matrix_sub >= 10) >= 5
count_matrix_filt <- count_matrix_sub[keep, ]

write.csv(count_matrix_sub, "data/processed/counts_matrix_matched.csv")
write.csv(count_matrix_filt, "data/processed/counts_matrix_filtered.csv")
write.csv(meta_matched, "data/processed/metadata_matched.csv", row.names = FALSE)

print(dim(count_matrix_filt))
print(table(meta_matched$disease_group))
