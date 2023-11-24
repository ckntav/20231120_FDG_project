setwd("/Users/chris/Desktop/20231120_FDG_project")

library(tidyverse)

#
tab_filepath <- file.path("input", "mapping", "MOUSE_10090_idmapping.dat")
raw_mapping <- read_tsv(tab_filepath, col_names = FALSE) %>% 
  set_names("UniProtKB-AC", "ID_type", "ID")

head(raw_mapping)
raw_mapping$ID_type %>% unique
raw_mapping$ID_type %>% unique %>% length

# exploration
raw_mapping$`UniProtKB-AC` %>% unique %>% length

raw_mapping %>% dplyr::filter(`UniProtKB-AC` %in% c("Q9CQV8-1", "Q9CQV8")) %>% print(n = 100)

raw_mapping %>% dplyr::filter(ID_type == "Gene_Name")

raw_mapping %>% dplyr::filter(ID_type == "Gene_Name") %>%
  group_by(ID) %>% tally %>% arrange(desc(n))

raw_mapping %>% dplyr::filter(ID_type == "Gene_Name") %>%
  group_by(ID) %>% tally %>% arrange(desc(n)) %>% 
  ggplot(aes(x = "mock", y = n)) +
  geom_jitter(width = 0.5) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic()

raw_mapping %>% dplyr::filter(ID_type == "Gene_Name") %>%
  group_by(ID) %>% tally %>% 
  pull(n) %>% summary

#
mapping_wide <- raw_mapping %>% 
  # dplyr::filter(!ID_type %in% c("GI", "Gene_Synonym")) %>% 
  dplyr::filter(ID_type %in% c("Gene_Name")) %>% 
  pivot_wider(names_from = ID_type,
              values_from = ID)

mapping_wide
