setwd("/Users/chris/Desktop/20231120_FDG_project")

library(tidyverse)
library(ComplexHeatmap)
library(eulerr)

#
male_model_rds_filename <- file.path("output", "rds", "male_model_df.rds")
female_model_rds_filename <- file.path("output", "rds", "female_model_df.rds")
mix_model_rds_filename <- file.path("output", "rds", "mix_model_df.rds")

#
male_model <- readRDS(male_model_rds_filename)
upreg_male <- male_model %>% dplyr::filter(pval <= 0.05, log2FC >= 1) %>% pull(UniProtKB_ID)
downreg_male <- male_model %>% dplyr::filter(pval <= 0.05, log2FC <= -1) %>% pull(UniProtKB_ID)

#
female_model <- readRDS(female_model_rds_filename)
upreg_female <- female_model %>% dplyr::filter(pval <= 0.05, log2FC >= 1) %>% pull(UniProtKB_ID)
downreg_female <- female_model %>% dplyr::filter(pval <= 0.05, log2FC <= -1) %>% pull(UniProtKB_ID)

#
mix_model <- readRDS(mix_model_rds_filename)
upreg_mix <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 1) %>% pull(UniProtKB_ID)
mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 1) %>% pull(UniProtKB_AC) %>% paste(collapse = " ")
downreg_mix <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC <= -1) %>% pull(UniProtKB_ID)

#
all_sign_proteins <- list("upreg_male" = upreg_male,
                          "downreg_male" = downreg_male,
                          "upreg_female" = upreg_female,
                          "downreg_female" = downreg_female,
                          "upreg_mix" = upreg_mix,
                          "downreg_mix" = downreg_mix)
sapply(all_sign_proteins, length)

#
ltm <- list_to_matrix(all_sign_proteins)
combMat <- make_comb_mat(ltm)

annot_top <- HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(combMat), 
                                                                   border = FALSE,
                                                                   gp = gpar(fill = "black"),
                                                                   height = unit(3, "cm")), 
                               "Size" = anno_text(comb_size(combMat),
                                                  rot = 0,
                                                  just = "center",
                                                  location = 0.25),
                               annotation_name_side = "left", annotation_name_rot = 0)
annot_right <- rowAnnotation("Set size" = anno_barplot(set_size(combMat), 
                                                       border = FALSE, 
                                                       gp = gpar(fill = "black"), 
                                                       width = unit(2, "cm")),
                             "Size" = anno_text(set_size(combMat)))

UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right,
      set_order = c("upreg_male", "upreg_female", "upreg_mix",
                    "downreg_male", "downreg_female", "downreg_mix"))

fit_all_DEGs <- euler(ltm, shape = "ellipse")
plot(fit_all_DEGs, labels = TRUE, legend = list(side = "bottom"),
     # quantities = list(type = c("counts", "percent")),
     quantities = list(type = c("counts")),
     fills = list(fill = c("#2B70AB", "#FFB027", "#3EA742", "#CD3301", "#9370DB", "grey")),
     edges = list(alpha = 0))

#
sign_proteins_nomix <- list("upreg_male" = upreg_male,
                          "downreg_male" = downreg_male,
                          "upreg_female" = upreg_female,
                          "downreg_female" = downreg_female)
sapply(sign_proteins_nomix, length)

#
ltm <- list_to_matrix(sign_proteins_nomix)
combMat <- make_comb_mat(ltm)

annot_top <- HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(combMat), 
                                                                   border = FALSE,
                                                                   gp = gpar(fill = "black"),
                                                                   height = unit(3, "cm")), 
                               "Size" = anno_text(comb_size(combMat),
                                                  rot = 0,
                                                  just = "center",
                                                  location = 0.25),
                               annotation_name_side = "left", annotation_name_rot = 0)
annot_right <- rowAnnotation("Set size" = anno_barplot(set_size(combMat), 
                                                       border = FALSE, 
                                                       gp = gpar(fill = "black"), 
                                                       width = unit(2, "cm")),
                             "Size" = anno_text(set_size(combMat)))

UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right,
      set_order = c("upreg_male", "upreg_female",
                    "downreg_male", "downreg_female"))

fit_nomix <- euler(ltm, shape = "ellipse")
plot(fit_nomix, labels = TRUE, legend = list(side = "bottom"),
     # quantities = list(type = c("counts", "percent")),
     quantities = list(type = c("counts")),
     fills = list(fill = c("#2B70AB", "#FFB027", "#3EA742", "#CD3301", "#9370DB", "grey")),
     edges = list(alpha = 0))