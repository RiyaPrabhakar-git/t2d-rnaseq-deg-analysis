# ==============================
# 01_metadata_processing.R
# Download and clean GEO metadata
# ==============================

source("scripts/00_setup.R")

gse <- getGEO("GSE164416", GSEMatrix = TRUE)
meta <- pData(gse[[1]])

write.csv(meta, "data/processed/metadata_raw.csv", row.names = FALSE)

meta_clean <- meta %>%
  mutate(
    disease_group = case_when(
      grepl("_T2D$", title, ignore.case = TRUE) ~ "T2D",
      grepl("_ND$", title, ignore.case = TRUE) ~ "ND",
      grepl("_IGT$", title, ignore.case = TRUE) ~ "IGT",
      grepl("_T3cD$", title, ignore.case = TRUE) ~ "T3cD",
      TRUE ~ NA_character_
    ),
    sample_id = sub(".*_(DP[0-9]+)_.*", "\\1", title)
  )

write.csv(meta_clean, "data/processed/metadata_all_groups.csv", row.names = FALSE)

meta_clean <- meta_clean %>%
  filter(disease_group %in% c("ND", "T2D"))

write.csv(meta_clean, "data/processed/metadata_clean.csv", row.names = FALSE)

print(table(meta_clean$disease_group))
print(head(meta_clean[, c("title", "sample_id", "disease_group")]))
