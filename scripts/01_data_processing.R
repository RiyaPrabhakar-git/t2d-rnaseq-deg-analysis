# Load libraries
library(tidyverse)
library(GEOquery)

# Load metadata
gse <- getGEO("GSE164416", GSEMatrix = TRUE)
meta <- pData(gse[[1]])

# Clean metadata
meta_clean <- meta %>%
  mutate(
    disease_group = case_when(
      grepl("_T2D$", title) ~ "T2D",
      grepl("_ND$", title) ~ "ND",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(disease_group %in% c("ND", "T2D"))

meta_clean$sample_id <- sub(".*_(DP[0-9]+)_.*", "\\1", meta_clean$title)
