setwd("/Users/chris/Downloads/Felipe_analysis")

library(tidyverse)

#
selected_colnames <- c("UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)",
                       "RefSeq", "GI", "PDB", "GO", "UniRef100", "UniRef90",
                       "UniRef50", "UniParc", "PIR", "NCBI-taxon", "MIM",
                       "UniGene", "PubMed", "EMBL", "EMBL-CDS", "Ensembl",
                       "Ensembl_TRS", "Ensembl_PRO", "Additional PubMed")

#
tab_filename <- "MOUSE_10090_idmapping_selected.tab"
raw_mapping <- read_tsv(tab_filename, col_names = FALSE) %>% 
  set_names(selected_colnames)
mapping_uniprot <- raw_mapping %>% dplyr::select("UniProtKB-ID")

"H2A1F_MOUSE" %in% mapping_uniprot$`UniProtKB-ID`



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
         uniprot = tmp %>% map(1) %>% unlist,
         protein_name_tmp = tmp %>% map(2) %>% unlist,
         protein_name = str_replace(pattern = "_MOUSE", replacement = "", protein_name_tmp)) %>% 
  dplyr::select(-tmp1, -tmp, -protein_name_tmp) %>% 
  dplyr::filter(!is.na(pval))


max(mix_model$mlog10pval)
min(mix_model$mlog10pval)

max(-log10(mix_model$pval)) + 5


#
library(EnhancedVolcano)

EnhancedVolcano(mix_model,
                lab = mix_model$protein_name,
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


sign_mixmodel <- mix_model %>% dplyr::filter(pval <= 0.05) %>% 
  mutate(direction = case_when(
    log2FC >= 0 ~ "upreg",
    log2FC <= 0 ~ "downreg"
  ))

sign_mixmodel %>% group_by(direction) %>% tally

sign_mixmodel %>% dplyr::filter(direction == "upreg") %>%
  arrange(pval) %>% 
  print(n = 101)

sign_mixmodel %>% dplyr::filter(direction == "downreg") %>%
  arrange(pval) %>% 
  print(n = 31)

#
candidate_protein <- c("EZRI")
candidate_uniprot <- c("P26040")
mix_model %>% dplyr::filter(protein_name %in% candidate_protein)
mix_model %>% dplyr::filter(uniprot %in% candidate_uniprot)
sign_mixmodel %>% dplyr::filter(uniprot %in% candidate_uniprot)

#
# Function to map UniProt IDs to protein names using UniProt REST API
map_uniprot_to_protein_names <- function(uniprot_ids) {
  base_url <- "https://www.ebi.ac.uk/proteins/api/proteins/"
  
  # Combine multiple UniProt IDs into a single string separated by commas
  id_string <- paste(uniprot_ids, collapse = ",")
  
  # Create the API request URL
  request_url <- paste0(base_url, id_string, "?format=json")
  
  # Send the API request
  response <- httr::GET(request_url)
  
  # Check if the request was successful
  if (httr::http_error(response)) {
    stop("Error in API request. Check your UniProt IDs and try again.")
  }
  
  # Parse the JSON response
  data <- httr::content(response, "parsed")
  
  # Extract protein names from the response
  protein_names <- sapply(data, function(entry) entry$name)
  
  # Return a data frame with UniProt IDs and corresponding protein names
  data.frame(UniProt_ID = uniprot_ids, Protein_Name = protein_names, stringsAsFactors = FALSE)
}

# Example usage
uniprot_ids <- c("P12345", "Q98765", "O54321")
result <- map_uniprot_to_protein_names(uniprot_ids)
print(result)
