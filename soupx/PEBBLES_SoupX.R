# This program will eliminate ambient RNA via SoupX

####################
## Load Libraries ##
####################
library(Seurat)
library(SoupX)
library(dplyr)

###########################
## Set working directory ##
###########################
setwd('/share/lasallelab/data/2020_Rett_Mouse_SingleCellRNAseq_Osman/03_soupx')

#############################
## Load Samples with SoupX ##
#############################

# Define a function to load, process, and rename the data
process_data <- function(sample_name, folder_path) {
  # Load data
  sc <- load10X(paste0(folder_path, sample_name, "/outs/"))
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  colnames(out) <- paste0(sample_name, "_", colnames(out))
  return(out)
}

# List of sample names
samples <- c(
  "24_PCB_WT",
  "25_VEHICLE_WT",
  "27_PCB_HET_2",
  "27_PCB_HET",
  "28_VEHICLE_HET",
  "29_VEHICLE_WT",
  "30_VEHICLE_WT_2",
  "30_VEHICLE_WT",
  "31_PCB_WT",
  "36_PCB_HET",
  "37_PCB_WT_2",
  "37_PCB_WT",
  "38_VEHICLE_HET",
  "39_PCB_HET",
  "40_VEHICLE_HET_2",
  "40_VEHICLE_HET"
)

# Folder path
folder_path <- "/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/"

###################
## Decontaminate ##
###################

# Process data for each sample
scs <- lapply(samples, function(sample_name) {
  process_data(sample_name, folder_path)
})


#Combine all count matricies into one matrix and create a seurat object
srat = do.call(cbind,scs)
PEBBLES_soupx = CreateSeuratObject(srat,
                          project = "PEBBLES_soupx",
                          min.cells = 10,
                          min.features = 200,
                          names.field = 2,
                          names.delim = "\\-")

# Add a new columns metadata
samples <- sub("_[^_]+$", "", rownames(PEBBLES_soupx@meta.data))
PEBBLES_soupx$meta.data$Samples <- samples
PEBBLES_soupx@meta.data$Genotype <- ifelse(grepl("WT", PEBBLES_soupx@meta.data$Samples), "WT",
                                           ifelse(grepl("HET", PEBBLES_soupx@meta.data$Samples), "HET", NA))
PEBBLES_soupx@meta.data$Treatment <- ifelse(grepl("PCB", PEBBLES_soupx@meta.data$Samples), "PCB",
                                           ifelse(grepl("VEHICLE", PEBBLES_soupx@meta.data$Samples), "VEHICLE", NA))

#################################
## Save the SoupX Seurat object##
#################################
save(PEBBLES_soupx,file="PEBBLES_soupx.RData")
