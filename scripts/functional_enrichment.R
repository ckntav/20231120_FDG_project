setwd("/Users/chris/Desktop/20231120_FDG_project")

library(tidyverse)
library(ComplexHeatmap)
library(eulerr)

#
mapping_ids <- readRDS(file = "output/rds/mapping_ids.rds")

#
mix_model_rds_filename <- file.path("output", "rds", "mix_model_df.rds")
mix_model <- readRDS(mix_model_rds_filename) %>% left_join(mapping_ids)

#
upreg_mix_df <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 0) %>% 
  dplyr::select(UniProtKB_ID, gene, `GeneID_(EntrezGene)`, Ensembl) %>% 
  print(n = 69)

upreg_mix_FC0 <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 0) %>% pull(UniProtKB_AC)
upreg_mix_FC1 <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 1) %>% pull(UniProtKB_AC)

# upreg_mix <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 1) %>% pull(`GeneID_(EntrezGene)`)

# mix_model %>% dplyr::filter(pval <= 0.05, log2FC >= 0) %>% pull(UniProtKB_AC) %>% paste(collapse = " ")

downreg_mix_FC0 <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC <= -0) %>% pull(UniProtKB_AC)
downreg_mix_FC1 <- mix_model %>% dplyr::filter(pval <= 0.05, log2FC <= -1) %>% pull(UniProtKB_AC)

#
# library(clusterProfiler)
# library(org.Mm.eg.db)
# library(AnnotationDbi)
# 
# GO_results_BP_upreg <- enrichGO(gene = upreg_mix, OrgDb = org.Mm.eg.db, ont = "BP")
# dotplot(GO_results_BP_upreg, showCategory = 50)
# 
# upreg_genes_entrezId <- clusterProfiler::bitr(geneID = upreg_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
# upreg_kk <- enrichKEGG(upreg_mix, organism = "mmu", minGSSize = 1)
# upreg_kk %>% as.data.frame

# library(ReactomePA)
# upregRPA <- enrichPathway(upreg_mix, organism = "mouse", minGSSize = 1, readable = TRUE)
# upregRPA %>% as.data.frame

#
# install.packages("gprofiler2")

library(gprofiler2)

#
# gostres_upreg <- gost(upreg_mix, organism = "mmusculus", correction_method = "gSCS")
# df_res_upreg <- gostres_upreg$result
# df_res_upreg %>%
#   dplyr::select(term_name, term_id, p_value, source, everything())
# p_upreg <- gostplot(gostres_upreg, capped = FALSE, interactive = TRUE)
# p_upreg
# 
# #
# gostres_downreg <- gost(downreg_mix, organism = "mmusculus", correction_method = "gSCS")
# df_res_downreg <- gostres_downreg$result
# df_res_downreg
# p_downreg <- gostplot(gostres_downreg, capped = FALSE, interactive = TRUE)
# p_downreg

#
diff_expressed_FC1 <- list()
diff_expressed_FC1$upregulated_proteins_FC1 <- upreg_mix_FC1
diff_expressed_FC1$downregulated_proteins_FC1 <- downreg_mix_FC1

gostres_diff_expressed_FC1 <- gost(diff_expressed_FC1, organism = "mmusculus", correction_method = "gSCS")
df_res_FC1 <- gostres_diff_expressed_FC1$result
p_FC1 <- gostplot(gostres_diff_expressed_FC1, capped = FALSE, interactive = TRUE)
p_FC1

write_tsv(df_res_FC1, file = "output/tsv/20231123_pathways_FC1.tsv")

#
diff_expressed_FC0 <- list()
diff_expressed_FC0$upregulated_proteins_FC0 <- upreg_mix_FC0
diff_expressed_FC0$downregulated_proteins_FC0 <- downreg_mix_FC0

gostres_diff_expressed_FC0 <- gost(diff_expressed_FC0, organism = "mmusculus", correction_method = "gSCS")
df_res_FC0 <- gostres_diff_expressed_FC0$result
p_FC0 <- gostplot(gostres_diff_expressed_FC0, capped = FALSE, interactive = TRUE)
p_FC0

write_tsv(df_res_FC0, file = "output/tsv/20231123_pathways_FC0.tsv")


