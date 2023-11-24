setwd("/Users/chris/Desktop/20231120_FDG_project")

library(tidyverse)
source("scripts/ckn_utils/ckn_utils_savePlot.R")

#
tsv_filepath <- file.path("input", "docs", "Steve Lacroix 2797 070823 volcano plot -tableau d'origine.txt")
raw_data <- read_tsv(tsv_filepath)

colnames(raw_data)
head(raw_data %>% as.data.frame)

#
female_model <- raw_data %>% 
  dplyr::select("T: Accession", "N: Log2 (KO #1 Female/WT #1 Female)", "N: -Log p-value (KO #1 Female/WT #1 Female)", "T: Description") %>% 
  set_names("accession_id", "log2FC_chr", "mlog10pval_chr", "description") %>% 
  mutate(log2FC = str_replace(pattern = ",", replacement = "\\.", string = log2FC_chr) %>% as.numeric,
         mlog10pval = str_replace(pattern = ",", replacement = "\\.", string = mlog10pval_chr) %>% as.numeric) %>% 
  dplyr::select(-log2FC_chr, -mlog10pval_chr) %>% 
  mutate(pval = 10^(-mlog10pval)) %>% 
  mutate(tmp1 = str_replace(pattern = "tr\\|", replacement = "", accession_id),
         tmp = str_split(string = tmp1, pattern = "\\|"),
         UniProtKB_AC = tmp %>% map(1) %>% unlist,
         UniProtKB_ID = tmp %>% map(2) %>% unlist,
         protein_name = str_replace(pattern = "_MOUSE", replacement = "", UniProtKB_ID)) %>% 
  dplyr::select(-tmp1, -tmp) %>% 
  dplyr::filter(!is.na(pval))

#
female_model_rds_filename <- file.path("output", "rds", "female_model_df.rds")
saveRDS(female_model, file = female_model_rds_filename)

#
max(female_model$mlog10pval)
min(female_model$mlog10pval)

max(-log10(female_model$pval)) + 5

#
library(EnhancedVolcano)

volcano_plot <-
EnhancedVolcano(female_model,
                lab = female_model$protein_name,
                x = "log2FC",
                y = "pval",
                FCcutoff = 1,
                pCutoff = 5e-2,
                pointSize = 0.5,
                labSize = 3,
                ylim = c(0, max(-log10(female_model$pval)) + 0.5),
                title = "IL-1β KO vs IL-1β WT",
                subtitle = "MALE KO / MALE WT",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = "bottom")
volcano_plot

sign_female_model <- female_model %>% dplyr::filter(pval <= 0.05) %>% 
  mutate(direction = case_when(
    log2FC >= 0 ~ "upreg",
    log2FC <= 0 ~ "downreg"
  ))

sign_female_model %>% group_by(direction) %>% tally

sign_female_model %>% dplyr::filter(direction == "upreg") %>%
  arrange(pval)

sign_female_model %>% dplyr::filter(direction == "downreg") %>%
  arrange(pval)

sign_female_model %>% dplyr::filter(direction == "upreg") %>% pull(UniProtKB_ID)
