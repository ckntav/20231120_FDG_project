setwd("/Users/chris/Desktop/20231120_FDG_project")

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
source("scripts/ckn_utils/ckn_utils_savePlot.R")

#
tsv_filepath <- file.path("input", "docs", "Steve Lacroix 2797 070823 volcano plot -tableau d'origine.txt")
raw_data <- read_tsv(tsv_filepath, locale = locale(decimal_mark = ","))

colnames(raw_data)
head(raw_data %>% as.data.frame)

#
ivalues_mix <- raw_data %>% 
  dplyr::select("T: Accession", "T: Description",
                "WT #1 (1) Area", "WT #1 (2) Area", "KO #1 (3) Area", "KO #1 (4) Area",
                "WT #2 (5) Area", "WT #2 (6) Area", "KO #2 (7) Area", "KO #2 (8) Area",
                "MIXC57 #3 (9) Area", "MIXC57 #3 (10) Area", "MIXKO #3 (11) Area", "MIXKO #3 (12) Area") %>% 
  set_names("accession_id", "description",
            "male_wt1", "male_wt2", "male_ko1", "male_ko2",
            "female_wt1", "female_wt2", "female_ko1", "female_ko2",
            "mix_wt1", "mix_wt2", "mix_ko1", "mix_ko2")

#
mat <- ivalues_mix %>% dplyr::select(-accession_id, -description) %>% as.matrix
colnames(mat) <- c("male_wt1", "male_wt2", "male_ko1", "male_ko2",
                   "female_wt1", "female_wt2", "female_ko1", "female_ko2",
                   "mix_wt1", "mix_wt2", "mix_ko1", "mix_ko2")

mat.pca <- prcomp(mat)
summary(mat.pca)

df <- cbind(
  data.frame("sample_id" = rownames(mat.pca$rotation)),
  mat.pca$rotation[, c(1,2)]
)

df %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text(label = df$sample_id, vjust = 2) +
  theme_minimal() +
  xlim(c(-0.305, -0.27)) +
  ylim(c(-0.6, 0.5))
