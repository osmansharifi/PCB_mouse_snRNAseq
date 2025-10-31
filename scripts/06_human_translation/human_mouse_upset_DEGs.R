##############################
## Upset plot of KEGG terms ##
##############################
# Author : Osman Sharifi

# Load libraries
library(UpSetR)
library(dplyr)
library(ggplot2)
library(glue)

#load data
human_degs <- read.csv('/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/Human_rett_PCBvs_noPCB_sig_DEGs.csv')
human_glut_list <- human_degs %>% filter(Cell_type == "Glutamatergic")
human_gaba_list <- human_degs %>% filter(Cell_type == "GABAergic")
human_nonn_list <- human_degs %>% filter(Cell_type == "Non-Neuronal")
mouse_degs <- read.csv('/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE.csv')
# Convert the SYMBOL column to uppercase
mouse_degs <- mouse_degs %>%
  mutate(SYMBOL = toupper(SYMBOL))
mouse_glut_list <- mouse_degs %>% filter(Cell_type == "Glutamatergic")
mouse_gaba_list <- mouse_degs %>% filter(Cell_type == "GABAergic")
mouse_nonn_list <- mouse_degs %>% filter(Cell_type == "Non-neuronal")

#######################
## Create upset plot ##
#######################
listInput <- list(Glutamatergic_human = human_glut_list$SYMBOL, 
                  Glutamatergic_mouse = mouse_glut_list$SYMBOL,
                  GABAergic_human = human_gaba_list$SYMBOL, 
                  GABAergic_mouse = mouse_gaba_list$SYMBOL,
                  Non_neuronal_human = human_nonn_list$SYMBOL,
                  Non_neuronal_mouse = mouse_nonn_list$SYMBOL)
# Find common SYMBOLs across all lists
common_symbols <- Reduce(intersect, listInput)

# View the result
common_symbols
# Find common SYMBOLs between GABAergic_human and GABAergic_mouse
common_gaba_symbols <- intersect(listInput$GABAergic_human, listInput$GABAergic_mouse)

newlist <- list(GABAergic_human = human_nonn_list$SYMBOL,  GABAergic_mouse = mouse_nonn_list$SYMBOL)
# View the result
common_gaba_symbols <- Reduce(intersect, newlist)
common_gaba_symbols

base_path <- "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/"
pdf(glue("{base_path}upset_human_mouse_DEGs.pdf"))
upset(fromList(listInput), sets = c('Glutamatergic_human','Glutamatergic_mouse', 'GABAergic_human', 'GABAergic_mouse', 'Non_neuronal_human', 'Non_neuronal_mouse'), keep.order = TRUE)
dev.off()

