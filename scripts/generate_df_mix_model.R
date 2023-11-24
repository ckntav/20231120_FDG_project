setwd("/Users/chris/Desktop/20231120_FDG_project")

library(tidyverse)
source("scripts/ckn_utils/ckn_utils_savePlot.R")

#
tsv_filepath <- file.path("input", "docs", "Steve Lacroix 2797 070823 volcano plot -tableau d'origine.txt")
raw_data <- read_tsv(tsv_filepath)

colnames(raw_data)
head(raw_data %>% as.data.frame)

#
mix_model <- raw_data %>% 
  dplyr::select("T: Accession", "N: Log2 (MIXKO/MIXWT)", "N: -Log p-value (MIXKO/MIXWT)", "T: Description") %>% 
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
  dplyr::filter(!is.na(pval)) %>% 
  mutate(tmp = str_split(pattern = " ", description) %>% map(-3) %>% unlist,
         gene = str_replace(string = tmp, pattern = "GN=", replacement = "")) %>% 
  dplyr::select(-tmp)

#
mix_model_rds_filename <- file.path("output", "rds", "mix_model_df.rds")
saveRDS(mix_model, file = mix_model_rds_filename)

#
max(mix_model$mlog10pval)
min(mix_model$mlog10pval)

max(-log10(mix_model$pval)) + 5

#
library(EnhancedVolcano)

volcano_plot <-
EnhancedVolcano(mix_model,
                lab = mix_model$gene,
                x = "log2FC",
                y = "pval",
                FCcutoff = 1,
                pCutoff = 5e-2,
                pointSize = 0.5,
                labSize = 3,
                ylim = c(0, max(-log10(mix_model$pval)) + 0.5),
                title = "IL-1β KO vs IL-1β WT",
                subtitle = "MIXKO/MIXWT",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = "bottom")

volcano_plot

sign_mix_model <- mix_model %>% dplyr::filter(pval <= 0.05) %>% 
  mutate(direction = case_when(
    log2FC >= 1 ~ "upreg",
    log2FC <= -1 ~ "downreg",
    TRUE ~ "cst"
  ))

sign_mix_model %>% group_by(direction) %>% tally

#
sign_mix_model %>% dplyr::filter(direction == "upreg") %>%
  arrange(pval) %>% 
  print(n = 69)

sign_mix_model %>% dplyr::filter(direction == "upreg") %>%
  arrange(pval) %>% pull(gene) %>% unique

sign_mix_model %>% dplyr::filter(direction == "upreg") %>%
  arrange(pval) %>% 
  dplyr::select(accession_id, description) %>% 
  print(n = 69)

#
sign_mix_model %>% dplyr::filter(direction == "downreg") %>%
  arrange(pval) %>% 
  print(n = 25)

sign_mix_model %>% dplyr::filter(direction == "downreg") %>%
  arrange(pval) %>% pull(gene) %>% unique

sign_mix_model %>% dplyr::filter(direction == "downreg") %>%
  arrange(pval) %>% 
  dplyr::select(accession_id, description) %>%
  print(n = 25)

#
upreg_mix_model <- sign_mix_model %>% dplyr::filter(direction == "upreg")
nrow(upreg_mix_model)
upreg_mix_model %>% pull(gene) %>% unique
write_tsv(upreg_mix_model, file = "output/tsv/20231123_upreg_mix_model.tsv")

downreg_mix_model <- sign_mix_model %>% dplyr::filter(direction == "downreg")
nrow(downreg_mix_model)
downreg_mix_model %>% pull(gene) %>% unique
write_tsv(downreg_mix_model, file = "output/tsv/20231123_downreg_mix_model.tsv")
