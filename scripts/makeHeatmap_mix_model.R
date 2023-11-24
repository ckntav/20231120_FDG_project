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
mix_model_rds_filename <- file.path("output", "rds", "mix_model_df.rds")
mix_model <- readRDS(mix_model_rds_filename)
proteins_of_interest <- mix_model %>% dplyr::filter(pval <= 0.05, abs(log2FC) >= 1) %>% pull(accession_id)
upreg_proteins <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 1) %>% pull(accession_id)
downreg_proteins <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC <= -1) %>% pull(accession_id)


#
ivalues_mix <- raw_data %>% 
  dplyr::select("T: Accession", "MIXC57 #3 (9) Area", "MIXC57 #3 (10) Area", "MIXKO #3 (11) Area", "MIXKO #3 (12) Area", "T: Description") %>% 
  set_names("accession_id", "mix_wt1_chr", "mix_wt2_chr", "mix_ko1_chr", "mix_ko2_chr", "description") %>% 
  mutate(mix_wt1 = str_replace(pattern = ",", replacement = "\\.", string = mix_wt1_chr) %>% as.numeric,
         mix_wt2 = str_replace(pattern = ",", replacement = "\\.", string = mix_wt2_chr) %>% as.numeric,
         mix_ko1 = str_replace(pattern = ",", replacement = "\\.", string = mix_ko1_chr) %>% as.numeric,
         mix_ko2 = str_replace(pattern = ",", replacement = "\\.", string = mix_ko2_chr) %>% as.numeric) %>% 
  mutate(tmp1 = str_replace(pattern = "tr\\|", replacement = "", accession_id),
         tmp = str_split(string = tmp1, pattern = "\\|"),
         UniProtKB_AC = tmp %>% map(1) %>% unlist,
         UniProtKB_ID = tmp %>% map(2) %>% unlist,
         protein_name = str_replace(pattern = "_MOUSE", replacement = "", UniProtKB_ID)) %>% 
  mutate(tmp = str_split(pattern = " ", description) %>% map(-3) %>% unlist,
         gene = str_replace(string = tmp, pattern = "GN=", replacement = "")) %>% 
  mutate(label = paste0(gene, "_", protein_name)) %>% 
  dplyr::select(accession_id, label, mix_wt1, mix_wt2, mix_ko1, mix_ko2) %>% 
  dplyr::filter(accession_id %in% proteins_of_interest) %>% 
  mutate(direction = case_when(
    accession_id %in% upreg_proteins ~ "upreg",
    accession_id %in% downreg_proteins ~ "downreg",
    TRUE ~ 'cst'
  )) %>% 
  arrange(desc(direction), mix_wt1, mix_wt2, mix_ko1, mix_ko2)

##### all upreg and downreg
mat_direction <- ivalues_mix %>% dplyr::select(direction) %>% as.matrix

#
annot_right <- rowAnnotation(direction = mat_direction[, 1],
                             col = list(direction = c("upreg" = "#EB5000",
                                                      "downreg" = "#ACC07E",
                                                      "cst" = "#E9C622")),
                             annotation_name_side = "top")
#
mat_mix <- ivalues_mix %>%
  dplyr::select(mix_wt1, mix_wt2, mix_ko1, mix_ko2) %>% 
  as.matrix
rownames(mat_mix) <- ivalues_mix$label

# icolors <- colorRamp2(c(13, 30), c("#DAA520", "white", "#79A9CD"))
# icolors <- colorRamp2(c(13, 30), c("#DAA520", "#79A9CD"))

Heatmap(mat_mix,
        name = "intensity",
        # col = icolors,
        # col = rev(magma(100)),
        # col = rev(mako(100)),
        col = rev(rocket(100)),
        right_annotation = annot_right,
        column_names_side = "top", column_names_rot = 45,
        row_names_gp = gpar(fontsize = 7),
        cluster_rows = FALSE,
        cluster_columns = FALSE)

# ##### only upreg
# mat_direction_upreg <- ivalues_mix %>% dplyr::filter(direction == "upreg") %>% dplyr::select(direction) %>% as.matrix
# 
# #
# annot_right_upreg <- rowAnnotation(direction = mat_direction_upreg[, 1],
#                              col = list(direction = c("upreg" = "#EB5000",
#                                                       "downreg" = "#ACC07E",
#                                                       "cst" = "#E9C622")))
#
# mat_mix_upreg <- ivalues_mix %>% dplyr::filter(direction == "upreg") %>% 
#   dplyr::select(mix_wt1, mix_wt2, mix_ko1, mix_ko2) %>% 
#   as.matrix
# rownames(mat_mix_upreg) <- ivalues_mix %>% dplyr::filter(direction == "upreg") %>% pull(label)
# 
# min(mat_mix_upreg)
# max(mat_mix_upreg)
# 
# icolors <- colorRamp2(c(13, 30), c("#DAA520", "white", "#79A9CD"))
# icolors <- colorRamp2(c(13, 29), c("white", "#DAA520"))
# 
# Heatmap(mat_mix_upreg,
#         name = "intensity",
#         # col = icolors,
#         col = viridis(100),
#         right_annotation = annot_right_upreg,
#         cluster_rows = FALSE,
#         cluster_columns = FALSE)
# 
