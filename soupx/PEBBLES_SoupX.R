# This program will eliminate ambient RNA via SoupX

####################
## Load Libraries ##
####################
library(Seurat)
library(SoupX)

###########################
## Set working directory ##
###########################
setwd('/share/lasallelab/data/2020_Rett_Mouse_SingleCellRNAseq_Osman/03_soupx')

#############################
## Load Samples with SoupX ##
#############################

# load all samples 
sc24_PCB_WT = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/24_PCB_WT/outs/')
sc24_PCB_WT = autoEstCont(sc24_PCB_WT)
out24_PCB_WT = adjustCounts(sc24_PCB_WT)
colnames(out24_PCB_WT) = paste0(nom,'_',colnames(out24_PCB_WT))
# Estimated global rho of 0.07
sc25_VEHICLE_WT = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/25_VEHICLE_WT/outs/')
sc25_VEHICLE_WT = autoEstCont(sc25_VEHICLE_WT)
out25_VEHICLE_WT = adjustCounts(sc25_VEHICLE_WT)
colnames(out25_VEHICLE_WT) = paste0(nom,'_',colnames(out25_VEHICLE_WT))
# Estimated global rho of 0.07
sc27_PCB_HET = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/27_PCB_HET/outs/')
sc27_PCB_HET = autoEstCont(sc27_PCB_HET)
out27_PCB_HET = adjustCounts(sc27_PCB_HET)
colnames(out27_PCB_HET) = paste0(nom,'_',colnames(out27_PCB_HET))
# Estimated global rho of 0.06
sc27_PCB_HET_2 = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/27_PCB_HET_2/outs/')
sc27_PCB_HET_2 = autoEstCont(sc27_PCB_HET_2)
out27_PCB_HET_2 = adjustCounts(sc27_PCB_HET_2)
colnames(out27_PCB_HET_2) = paste0(nom,'_',colnames(out27_PCB_HET_2))
# Estimated global rho of 0.04
sc28_VEHICLE_HET = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/28_VEHICLE_HET/outs/')
sc28_VEHICLE_HET = autoEstCont(sc28_VEHICLE_HET)
out28_VEHICLE_HET = adjustCounts(sc28_VEHICLE_HET)
colnames(out28_VEHICLE_HET) = paste0(nom,'_',colnames(out28_VEHICLE_HET))
# Estimated global rho of 0.07
sc29_VEHICLE_WT = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/29_VEHICLE_WT/outs/')
sc29_VEHICLE_WT = autoEstCont(sc29_VEHICLE_WT)
out28_VEHICLE_HET = adjustCounts(sc29_VEHICLE_WT)
colnames(out28_VEHICLE_HET) = paste0(nom,'_',colnames(out28_VEHICLE_HET))
# Estimated global rho of 0.06
sc30_VEHICLE_WT = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/30_VEHICLE_WT/outs/')
sc30_VEHICLE_WT = autoEstCont(sc30_VEHICLE_WT)
out30_VEHICLE_WT = adjustCounts(sc30_VEHICLE_WT)
colnames(out30_VEHICLE_WT) = paste0(nom,'_',colnames(out30_VEHICLE_WT))
# Estimated global rho of 0.04
sc30_VEHICLE_WT_2 = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/30_VEHICLE_WT_2/outs/')
sc30_VEHICLE_WT_2 = autoEstCont(sc30_VEHICLE_WT_2)
out30_VEHICLE_WT_2 = adjustCounts(sc30_VEHICLE_WT_2)
colnames(out30_VEHICLE_WT_2) = paste0(nom,'_',colnames(out30_VEHICLE_WT_2))
# Estimated global rho of 0.07
sc31_PCB_WT = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/31_PCB_WT/outs/')
sc31_PCB_WT = autoEstCont(sc31_PCB_WT)
out31_PCB_WT = adjustCounts(sc31_PCB_WT)
colnames(out31_PCB_WT) = paste0(nom,'_',colnames(out31_PCB_WT))
# Estimated global rho of 0.07
sc36_PCB_HET = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/36_PCB_HET/outs/')
sc36_PCB_HET = autoEstCont(sc36_PCB_HET)
out36_PCB_HET = adjustCounts(sc36_PCB_HET)
colnames(out36_PCB_HET) = paste0(nom,'_',colnames(out36_PCB_HET))
# Estimated global rho of 0.07
sc37_PCB_WT = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/37_PCB_WT/outs/')
sc37_PCB_WT = autoEstCont(sc37_PCB_WT)
out37_PCB_WT = adjustCounts(sc37_PCB_WT)
colnames(out37_PCB_WT) = paste0(nom,'_',colnames(out37_PCB_WT))
# Estimated global rho of 0.07
sc37_PCB_WT_2 = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/37_PCB_WT_2/outs/')
sc37_PCB_WT_2 = autoEstCont(sc37_PCB_WT_2)
out37_PCB_WT_2 = adjustCounts(sc37_PCB_WT_2)
colnames(out37_PCB_WT_2) = paste0(nom,'_',colnames(out37_PCB_WT_2))
# Estimated global rho of 0.06
sc38_VEHICLE_HET = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/38_VEHICLE_HET/outs/')
sc38_VEHICLE_HET = autoEstCont(sc38_VEHICLE_HET)
out38_VEHICLE_HET = adjustCounts(sc38_VEHICLE_HET)
colnames(out38_VEHICLE_HET) = paste0(nom,'_',colnames(out38_VEHICLE_HET))
# Estimated global rho of 0.04
sc39_PCB_HET = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/39_PCB_HET/outs/')
sc39_PCB_HET = autoEstCont(sc39_PCB_HET)
out39_PCB_HET = adjustCounts(sc39_PCB_HET)
colnames(out39_PCB_HET) = paste0(nom,'_',colnames(out39_PCB_HET))
# Estimated global rho of 0.04
sc40_VEHICLE_HET = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/40_VEHICLE_HET/outs/')
sc40_VEHICLE_HET = autoEstCont(sc40_VEHICLE_HET)
out40_VEHICLE_HET = adjustCounts(sc40_VEHICLE_HET)
colnames(out40_VEHICLE_HET) = paste0(nom,'_',colnames(out40_VEHICLE_HET))
# Estimated global rho of 0.08
sc40_VEHICLE_HET_2 = load10X('/share/lasallelab/Osman/2021_PEBBLES_Cortex/2021_mouse_PCB_raw_reads/01-Cellranger/40_VEHICLE_HET_2/outs/')
sc40_VEHICLE_HET_2 = autoEstCont(sc40_VEHICLE_HET_2)
out40_VEHICLE_HET_2 = adjustCounts(sc40_VEHICLE_HET_2)
colnames(out40_VEHICLE_HET_2) = paste0(nom,'_',colnames(out40_VEHICLE_HET_2))
# Estimated global rho of 0.05

###################
## Decontaminate ##
###################

scs = list(out24_PCB_WT, 
           out25_VEHICLE_WT, 
           out27_PCB_HET_2, 
           out27_PCB_HET,
           out28_VEHICLE_HET,
           out29_VEHICLE_WT,
           out30_VEHICLE_WT_2,
           out30_VEHICLE_WT,
           out31_PCB_WT,
           out36_PCB_HET,
           out37_PCB_WT_2,
           out37_PCB_WT,
           out38_VEHICLE_HET,
           out39_PCB_HET,
           out40_VEHICLE_HET_2,
           out40_VEHICLE_HET)

srat = list()
for(nom in names(scs)){
  #Clean channel named 'nom'
  tmp = adjustCounts(scs[[nom]])
  #Add experiment name to cell barcodes to make them unique
  colnames(tmp) = paste0(nom,'_',colnames(tmp))
  #Store the result
  srat[[nom]] = tmp
}
#Combine all count matricies into one matrix
srat = do.call(cbind,srat)
PEBBLES_soupx = CreateSeuratObject(srat,
                          project = "PEBBLES_soupx",
                          min.cells = 10,
                          min.features = 200,
                          names.field = 2,
                          names.delim = "\\-")

#################################
## Save the SoupX Seurat object##
#################################
save(PEBBLES_soupx,file="PEBBLES_soupx.RData")

