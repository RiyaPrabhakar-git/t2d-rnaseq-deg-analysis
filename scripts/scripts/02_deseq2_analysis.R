library(DESeq2)

counts <- read.delim("data/raw/GSE164416_DP_htseq_counts.txt")

# Prepare matrix
gene_col <- colnames(counts)[1]
count_matrix <- counts[, -1]
rownames(count_matrix) <- counts[[gene_col]]

# Matching
common_samples <- intersect(colnames(count_matrix), meta_clean$sample_id)

count_matrix_sub <- count_matrix[, common_samples]

meta_matched <- meta_clean[match(common_samples, meta_clean$sample_id), ]

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix_sub),
  colData = meta_matched,
  design = ~ disease_group
)

dds <- DESeq(dds)
