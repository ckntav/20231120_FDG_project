setwd("/Users/chris/Desktop/20231120_FDG_project")

library(tidyverse)

#
selected_colnames <- c("UniProtKB_AC", "UniProtKB_ID", "GeneID_(EntrezGene)",
                       "RefSeq", "GI", "PDB", "GO", "UniRef100", "UniRef90",
                       "UniRef50", "UniParc", "PIR", "NCBI_taxon", "MIM",
                       "UniGene", "PubMed", "EMBL", "EMBL_CDS", "Ensembl",
                       "Ensembl_TRS", "Ensembl_PRO", "Additional_PubMed")


tab_filename <- "input/mapping/MOUSE_10090_idmapping_selected.tab"
raw_mapping <- read_tsv(tab_filename, col_names = FALSE) %>% 
  set_names(selected_colnames)

mapping_ids <- raw_mapping %>% 
  dplyr::select("UniProtKB_AC", "UniProtKB_ID", "GeneID_(EntrezGene)", "Ensembl")

mapping_ids %>% group_by(Ensembl) %>% tally %>% arrange(desc(n))

saveRDS(mapping_ids, file = "output/rds/mapping_ids.rds")
