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
