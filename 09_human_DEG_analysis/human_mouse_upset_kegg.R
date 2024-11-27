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
file_paths <- list(
  Glut_human = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/Glutamatergic_enrichr.xlsx",
  GABA_human = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/GABAergic_enrichr.xlsx",
  NonN_human = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/Non-Neuronal_enrichr.xlsx",
  Glut_mouse = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/Glutamatergic_enrichr.xlsx",
  GABA_mouse = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/GABAergic_enrichr.xlsx",
  NonN_mouse = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/Non-neuronal_enrichr.xlsx"
)

output_files <- list(
  Glut_human = "Glut_human_filtered.xlsx",
  GABA_human = "GABA_human_filtered.xlsx",
  NonN_human = "NonN_human_filtered.xlsx",
  Glut_mouse = "Glut_mouse_filtered.xlsx",
  GABA_mouse = "GABA_mouse_filtered.xlsx",
  NonN_mouse = "NonN_mouse_filtered.xlsx"
)

# Function to load, filter, and save
process_file <- function(file_path, output_file, sheet = 4) {
  data <- read.xlsx(file_path, sheet = sheet) %>% 
    filter(P.value <= 0.05)
  write.xlsx(data, output_file)
  return(data)
}

# Apply to all datasets
filtered_data <- lapply(names(file_paths), function(name) {
  process_file(file_paths[[name]], output_files[[name]])
})

# Assign results to named list
names(filtered_data) <- names(file_paths)


#######################
## Create upset plot ##
#######################
listInput <- list(Glutamatergic_human = filtered_data$Glut_human$Term, 
                  Glutamatergic_mouse = filtered_data$Glut_mouse$Term,
                  GABAergic_human = filtered_data$GABA_human$Term, 
                  GABAergic_mouse = filtered_data$GABA_mouse$Term,
                  Non_neuronal_human = filtered_data$NonN_human$Term,
                  Non_neuronal_mouse = filtered_data$NonN_mouse$Term)
base_path <- "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/"
pdf(glue("{base_path}upset_human_mouse_kegg.pdf"))
upset(fromList(listInput), sets = c('Glutamatergic_human','Glutamatergic_mouse', 'GABAergic_human', 'GABAergic_mouse', 'Non_neuronal_human', 'Non_neuronal_mouse'), keep.order = TRUE)
dev.off()

