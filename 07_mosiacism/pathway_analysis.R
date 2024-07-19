#load all libraries
library(dplyr)
library(purrr)
library(openxlsx)
library(glue)
library(enrichR)

#load file containing all DEGs
all_DEGs <- read.csv('all_experiments_sig_DEGs.csv')

# Get unique combinations of DEG_experiment and Cell_type
combinations <- unique(select(all_DEGs, DEG_experiment, Cell_type))

# Loop through each combination
for (i in 1:nrow(combinations)) {
  DEG_experiment <- combinations$DEG_experiment[i]
  Cell_type <- combinations$Cell_type[i]
  
  # Subset data for the current combination
  subset_data <- filter(all_DEGs, DEG_experiment == DEG_experiment & Cell_type == Cell_type)
  
  # Ensure X column uniqueness
  subset_data <- subset_data %>%
    mutate(X = make.unique(X))  # Make X column values unique
  
  # Perform data manipulation steps
  DEGs <- subset_data %>%
    tibble::rownames_to_column(var = "RowID") %>%
    mutate(FC = case_when(logFC > 0 ~ 2^logFC,
                          logFC < 0 ~ -1/(2^logFC))) %>%
    select(RowID, SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %>%
    filter(P.Value < 0.05)
  
  # Write DEGs.xlsx
  write.xlsx(DEGs, file = glue("DEGs_{DEG_experiment}_{Cell_type}.xlsx"))
  
  # Filter sig DEGs and write sig_DEGs.xlsx
  sig_DEGs <- DEGs %>%
    filter(P.Value < 0.05) %>%
    select(-logFC)  # Assuming logFC is not required in sig_DEGs
  
  write.xlsx(sig_DEGs, file = glue("sig_DEGs_{DEG_experiment}_{Cell_type}.xlsx"))
  
  # Perform enrichment analysis
  enrichment_results <- DEGs %>%
    select(SYMBOL) %>%
    flatten() %>%
    enrichr(c("GO_Biological_Process_2023",
              "GO_Molecular_Function_2023",
              "GO_Cellular_Component_2023",
              "KEGG_2019_Mouse",
              "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
    set_names(names(.) %>% str_trunc(31, ellipsis = ""))
  
  # Write enrichment results to enrichr.xlsx
  write.xlsx(enrichment_results, file = glue("enrichr_{DEG_experiment}_{Cell_type}.xlsx"))
}
