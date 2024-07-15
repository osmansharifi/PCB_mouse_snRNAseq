######################
## Overlap analysis ##
######################
library(magrittr)
library(VennDiagram)
library(grDevices)
library(dplyr)
library(Seurat)
library(glue)
library(scCustomize)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(GeneOverlap)
library(enrichR)
library(readr)
################
## load files ##
################
# Set the working directory to the parent folder containing the subdirectories
setwd("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/")

# Get the list of directories
exp1 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/1_AllcellsVsAllcells_from_MUTPCB_WTPCB/sig_DEGs_1_AllcellsVsAllcells_from_MUTPCB_WTPCB.csv", header = TRUE)
exp2 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE.csv", header = TRUE)
exp3 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/3_WTcellsVsWTcells_from_WTPCB_WTVEHICLE/sig_WTcellsVsWTcells_from_WTPCB_WTVEHICLE_DEGs.csv", header = TRUE)
exp4 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/4_WTcellsVsWTcells_from_MUTPCB_WTPCB/sig_DEGs_4_WTcellsVsWTcells_from_MUTPCB_WTPCB.csv", header = TRUE)
exp5 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/5_WTcellsVsWTcells_from_MUTPCB_WTVEHICLE/sig_DEGs_5_WTcellsVsWTcells_from_MUTPCB_WTVEHICLE.csv", header = TRUE)
exp6 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/6_AllcellsVsAllcells_from_MUTVEHICLE_WTVEHICLE/sig_DEGs_6_AllcellsVsAllcells_from_MUTVEHICLE_WTVEHICLE.csv", header = TRUE)
exp7 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/7_WTcellsVsWTcells_from_MUTVEHICLE_WTVEHICLE/sig_DEGs_7_WTcellsVsWTcells_from_MUTVEHICLE_WTVEHICLE.csv", header = TRUE)
exp8 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/8_WTcellsVsWTcells_from_MUTVEHICLE_WTPCB/sig_DEGs_8_WTcellsVsWTcells_from_MUTVEHICLE_WTPCB.csv", header = TRUE)
exp9 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/9_WTcellsVsWTcells_from_HETPCB_HETVEHICLE/sig_WTcellsVsWTcells_from_HETPCB_HETVEHICLE_DEGs.csv", header = TRUE)
exp10 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/10_MUTcellsVsMUTcells_from_HETPCB_HETVEHICLE/sig_MUTcellsVsMUTcells_from_HETPCB_HETVEHICLE_DEGs.csv", header = TRUE)
exp11 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/11_WTcellsVsMUTcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_11_WTcellsVsMUTcells_from_MUTPCB_MUTVEHICLE.csv", header = TRUE)
exp12 <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/12_MUTcellsVsWTcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_12_MUTcellsVsWTcells_from_MUTPCB_MUTVEHICLE.csv", header = TRUE)

# List of file paths
file_paths <- c(
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/1_AllcellsVsAllcells_from_MUTPCB_WTPCB/sig_DEGs_1_AllcellsVsAllcells_from_MUTPCB_WTPCB.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/3_WTcellsVsWTcells_from_WTPCB_WTVEHICLE/sig_WTcellsVsWTcells_from_WTPCB_WTVEHICLE_DEGs.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/4_WTcellsVsWTcells_from_MUTPCB_WTPCB/sig_DEGs_4_WTcellsVsWTcells_from_MUTPCB_WTPCB.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/5_WTcellsVsWTcells_from_MUTPCB_WTVEHICLE/sig_DEGs_5_WTcellsVsWTcells_from_MUTPCB_WTVEHICLE.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/6_AllcellsVsAllcells_from_MUTVEHICLE_WTVEHICLE/sig_DEGs_6_AllcellsVsAllcells_from_MUTVEHICLE_WTVEHICLE.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/7_WTcellsVsWTcells_from_MUTVEHICLE_WTVEHICLE/sig_DEGs_7_WTcellsVsWTcells_from_MUTVEHICLE_WTVEHICLE.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/8_WTcellsVsWTcells_from_MUTVEHICLE_WTPCB/sig_DEGs_8_WTcellsVsWTcells_from_MUTVEHICLE_WTPCB.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/9_WTcellsVsWTcells_from_HETPCB_HETVEHICLE/sig_WTcellsVsWTcells_from_HETPCB_HETVEHICLE_DEGs.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/10_MUTcellsVsMUTcells_from_HETPCB_HETVEHICLE/sig_MUTcellsVsMUTcells_from_HETPCB_HETVEHICLE_DEGs.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/11_WTcellsVsMUTcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_11_WTcellsVsMUTcells_from_MUTPCB_MUTVEHICLE.csv",
  "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/12_MUTcellsVsWTcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_12_MUTcellsVsWTcells_from_MUTPCB_MUTVEHICLE.csv"
)

# Read each file into a list of data frames
exp_list <- lapply(file_paths, read.csv, header = TRUE)

# Extract the filenames without extension
file_names <- tools::file_path_sans_ext(basename(file_paths))

# Add DEG_experiment column to each data frame
for (i in seq_along(exp_list)) {
  exp_list[[i]]$DEG_experiment <- file_names[i]
}

combined_df <- do.call(rbind, exp_list)

# Write the concatenated data frame to a new CSV file
write_csv(concatenated_df, "concatenated_data.csv")

###########################################
## Venn diagram of the overlapping genes ##
###########################################

# Extract significant DEGs
sig_genes_Glutamatergic <- rownames(deg_results$Glutamatergic)[deg_results$Glutamatergic$adj.P.Val < 0.05]
sig_genes_GABAergic <- rownames(deg_results$GABAergic)[deg_results$GABAergic$adj.P.Val < 0.05]
sig_genes_Non_neuronal <- rownames(deg_results$`Non-neuronal`)[deg_results$`Non-neuronal`$adj.P.Val < 0.05]
intersection_all2 <- intersect(sig_genes_Non_neuronal,sig_genes_Glutamatergic)
intersection_all3 <- intersect(intersection_all2,sig_genes_GABAergic)
intersection_all4 <- intersect(sig_genes_Glutamatergic,sig_genes_GABAergic)
# Create a Venn diagram
pdf(glue("{directory_path}/venn_sig_WTcellsVsWTcells_from_WTPCB_WTVEHICLE_DEGs.pdf"))
temp <- venn.diagram(
  x = list(
    Glutamatergic = sig_genes_Glutamatergic,
    GABAergic = sig_genes_GABAergic,
    Non_neuronal = sig_genes_Non_neuronal
  ),
  category.names = c("Glutamatergic", "GABAergic", "Non-neuronal"),
  main = 'sig_WTcellsVsWTcells_from_HETPCB_HETVEHICLE_DEGs ',
  #filename = glue("{base_path}/broad_group_analysis/venn_glutamatergic.pdf"),
  filename = NULL,
  col = c('#E6B8BFFF', '#CC7A88FF', '#990F26FF'), 
  fill = c('#E6B8BFFF', '#CC7A88FF', '#990F26FF'),
  cat.cex = 1.2,
  cat.fontface = "bold",
  euler.d = TRUE,
  disable.logging = TRUE,
  hyper.test = TRUE
)
grid.draw(temp)
dev.off()


# Make all vectors the same length (pad with NA)
max_length <- max(length(sig_genes_Glutamatergic), length(sig_genes_GABAergic), length(sig_genes_Non_neuronal))
sig_genes_Glutamatergic <- c(sig_genes_Glutamatergic, rep(NA, max_length - length(sig_genes_Glutamatergic)))
sig_genes_GABAergic <- c(sig_genes_GABAergic, rep(NA, max_length - length(sig_genes_GABAergic)))
sig_genes_Non_neuronal <- c(sig_genes_Non_neuronal, rep(NA, max_length - length(sig_genes_Non_neuronal)))

# Combine into a dataframe
combined_df <- data.frame(sig_genes_Glutamatergic, sig_genes_GABAergic, sig_genes_Non_neuronal)

go.obj <- newGeneOverlap(combined_df$sig_genes_Glutamatergic,
                         combined_df$sig_genes_GABAergic,
                         combined_df$sig_genes_Non_neuronal,
                         genome.size = 22000)
go.obj <- testGeneOverlap(go.obj)
getPval(go.obj)
getOddsRatio(go.obj)
getJaccard(go.obj)
getContbl(go.obj)
show(go.obj)
writeLines(capture.output(show(go.obj)), glue('{directory_path}/geneoverlap_statistics.txt'))
