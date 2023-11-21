setwd("/Users/chris/Downloads/Felipe_analysis")

library(tidyverse)

#
tab_filename <- "MOUSE_10090_idmapping.dat"
raw_mapping <- read_tsv(tab_filename, col_names = FALSE) %>% 
  set_names("UniProtKB-AC", "ID_type", "ID")

# how many uniprot_ids?
raw_mapping$`UniProtKB-AC` %>% unique %>% length
raw_mapping %>% dplyr::filter(`UniProtKB-AC` %in% c("Q9CQV8-1", "Q9CQV8")) %>% print(n = 100)

raw_mapping %>% dplyr::filter(ID_type == "Gene_Name")