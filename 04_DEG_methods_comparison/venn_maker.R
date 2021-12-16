library(ggplot2)
library(ggVennDiagram)
library(Rcpp)
library(glue)

# Paths
Osman_limma_DEG_dir <- "~/Documents/GitHub/PEBBLES_mouse_snRNAseq/04_DEG_methods_comparison/P150_Female_Osman_Limma_DEGs/"
Kari_limma_DEG_dir <- "~/Documents/GitHub/PEBBLES_mouse_snRNAseq/04_DEG_methods_comparison/P150_Female_Kari_Limma_DEGs/"

venn_dir <- "~/Documents/GitHub/PEBBLES_mouse_snRNAseq/04_DEG_methods_comparison"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_F_P150_CORT"
subtitle_info <- "Mice, Female, P150, Cortex"

################################################################################
# Venn Diagram for Differentially Expressed Genes Per Analysis
# By Viktoria Haghani

for (cell_type in cell_types){
  # Read in data from Osman Limma analysis
  assign(paste0(cell_type, "_Osman_Limma_DEG"), read.csv(file = glue(Osman_limma_DEG_dir, cell_type, "_", metadata_info, "_Limma_DEG.csv")))
  # Read in data from Kari Limma analysis
  assign(paste0(cell_type, "_Kari_Limma_DEG"), read.csv(file = glue(Kari_limma_DEG_dir, cell_type, "_", metadata_info, "_Limma_DEG.csv")))
}

  # List of genes differentially expressed per cluster for Osman Limma
L2_3_IT_Osman_Limma_gene_list <- L2_3_IT_Osman_Limma_DEG$X
L6_Osman_Limma_gene_list <- L6_Osman_Limma_DEG$X
Sst_Osman_Limma_gene_list <- Sst_Osman_Limma_DEG$X
L5_Osman_Limma_gene_list <- L5_Osman_Limma_DEG$X
L4_Osman_Limma_gene_list <- L4_Osman_Limma_DEG$X
Pvalb_Osman_Limma_gene_list <- Pvalb_Osman_Limma_DEG$X
Sncg_Osman_Limma_gene_list <- Sncg_Osman_Limma_DEG$X
Non_neuronal_Osman_Limma_gene_list <- Non_neuronal_Osman_Limma_DEG$X
Oligo_Osman_Limma_gene_list <- Oligo_Osman_Limma_DEG$X
Vip_Osman_Limma_gene_list <- Vip_Osman_Limma_DEG$X
Lamp5_Osman_Limma_gene_list <- Lamp5_Osman_Limma_DEG$X
Astro_Osman_Limma_gene_list <- Astro_Osman_Limma_DEG$X
Peri_Osman_Limma_gene_list <- Peri_Osman_Limma_DEG$X
Endo_Osman_Limma_gene_list <- Endo_Osman_Limma_DEG$X
unique_Osman_Limma_genes <- unique(c(L2_3_IT_Osman_Limma_DEG$X,
                               L6_Osman_Limma_DEG$X,
                               Sst_Osman_Limma_DEG$X,
                               L5_Osman_Limma_DEG$X, 
                               L4_Osman_Limma_DEG$X,
                               Pvalb_Osman_Limma_DEG$X,
                               Sncg_Osman_Limma_DEG$X,
                               Non_neuronal_Osman_Limma_DEG$X,
                               Oligo_Osman_Limma_DEG$X,
                               Vip_Osman_Limma_DEG$X,
                               Lamp5_Osman_Limma_DEG$X,
                               Astro_Osman_Limma_DEG$X,
                               Peri_Osman_Limma_DEG$X,
                               Endo_Osman_Limma_DEG$X))

# List of genes differentially expressed per cluster for Kari Limma
L2_3_IT_Kari_Limma_gene_list <- L2_3_IT_Kari_Limma_DEG$SYMBOL
L6_Kari_Limma_gene_list <- L6_Kari_Limma_DEG$SYMBOL
Sst_Kari_Limma_gene_list <- Sst_Kari_Limma_DEG$SYMBOL
L5_Kari_Limma_gene_list <- L5_Kari_Limma_DEG$SYMBOL
L4_Kari_Limma_gene_list <- L4_Kari_Limma_DEG$SYMBOL
Pvalb_Kari_Limma_gene_list <- Pvalb_Kari_Limma_DEG$SYMBOL
Sncg_Kari_Limma_gene_list <- Sncg_Kari_Limma_DEG$SYMBOL
Non_neuronal_Kari_Limma_gene_list <- Non_neuronal_Kari_Limma_DEG$SYMBOL
Oligo_Kari_Limma_gene_list <- Oligo_Kari_Limma_DEG$SYMBOL
Vip_Kari_Limma_gene_list <- Vip_Kari_Limma_DEG$SYMBOL
Lamp5_Kari_Limma_gene_list <- Lamp5_Kari_Limma_DEG$SYMBOL
Astro_Kari_Limma_gene_list <- Astro_Kari_Limma_DEG$SYMBOL
Peri_Kari_Limma_gene_list <- Peri_Kari_Limma_DEG$SYMBOL
Endo_Kari_Limma_gene_list <- Endo_Kari_Limma_DEG$SYMBOL
unique_Kari_Limma_genes <- unique(c(L2_3_IT_Kari_Limma_DEG$SYMBOL,
                                     L6_Kari_Limma_DEG$SYMBOL,
                                     Sst_Kari_Limma_DEG$SYMBOL,
                                     L5_Kari_Limma_DEG$SYMBOL, 
                                     L4_Kari_Limma_DEG$SYMBOL,
                                     Pvalb_Kari_Limma_DEG$SYMBOL,
                                     Sncg_Kari_Limma_DEG$SYMBOL,
                                     Non_neuronal_Kari_Limma_DEG$SYMBOL,
                                     Oligo_Kari_Limma_DEG$SYMBOL,
                                     Vip_Kari_Limma_DEG$SYMBOL,
                                     Lamp5_Kari_Limma_DEG$SYMBOL,
                                     Astro_Kari_Limma_DEG$SYMBOL,
                                     Peri_Kari_Limma_DEG$SYMBOL,
                                     Endo_Kari_Limma_DEG$SYMBOL))
# Venn Diagram for Limma vs. DESeq2 vs. EdgeR per cluster
L2_3_IT_venn_list <- list(L2_3_IT_Osman_Limma_gene_list, L2_3_IT_Kari_Limma_gene_list)
L2_3_IT_venn <-ggVennDiagram(L2_3_IT_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma_Osman", "Limma_Kari")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L2_3_IT", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L2_3_IT_M_MUT_and_WT_F_P150_CORT_venn.pdf", device = "pdf", path = venn_dir)

L6_venn_list <- list(L6_Osman_Limma_gene_list, L6_Kari_Limma_gene_list)
L6_venn <- ggVennDiagram(L6_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma_Osman", "Limma_Kari")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L6", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L6_M_MUT_and_WT_F_P150_CORT_venn.pdf", device = "pdf", path = venn_dir)

Sst_venn_list <- list(Sst_Limma_gene_list, Sst_EdgeR_gene_list, Sst_DESeq2_gene_list)
Sst_venn <- ggVennDiagram(Sst_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Sst", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Sst_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

L5_venn_list <- list(L5_Limma_gene_list, L5_EdgeR_gene_list, L5_DESeq2_gene_list)
L5_venn <- ggVennDiagram(L5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L5", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L5_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

L4_venn_list <- list(L4_Limma_gene_list, L4_EdgeR_gene_list, L4_DESeq2_gene_list)
L4_venn <- ggVennDiagram(L4_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L4", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L4_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Pvalb_venn_list <- list(Pvalb_Limma_gene_list, Pvalb_EdgeR_gene_list, Pvalb_DESeq2_gene_list)
Pvalb_venn <- ggVennDiagram(Pvalb_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Pvalb", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Pvalb_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Sncg_venn_list <- list(Sncg_Limma_gene_list, Sncg_EdgeR_gene_list, Sncg_DESeq2_gene_list)
Sncg_venn <- ggVennDiagram(Sncg_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Sncg", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Sncg_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Non_neuronal_venn_list <- list(Non_neuronal_Limma_gene_list, Non_neuronal_EdgeR_gene_list, Non_neuronal_DESeq2_gene_list)
Non_neuronal_venn <- ggVennDiagram(Non_neuronal_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Non_neuronal", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Non_neuronal_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Oligo_venn_list <- list(Oligo_Limma_gene_list, Oligo_EdgeR_gene_list)
Oligo_venn <- ggVennDiagram(Oligo_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Oligo", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Oligo_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Vip_venn_list <- list(Vip_Limma_gene_list, Vip_EdgeR_gene_list, Vip_DESeq2_gene_list)
Vip_venn <- ggVennDiagram(Vip_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Vip", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Vip_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Lamp5_venn_list <- list(Lamp5_Limma_gene_list, Lamp5_EdgeR_gene_list, Lamp5_DESeq2_gene_list)
Lamp5_venn <- ggVennDiagram(Lamp5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Lamp5", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Lamp5_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Astro_venn_list <- list(Astro_Limma_gene_list, Astro_EdgeR_gene_list, Astro_DESeq2_gene_list)
Astro_venn <- ggVennDiagram(Astro_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Astro", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Astro_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Peri_venn_list <- list(Peri_Limma_gene_list, Peri_EdgeR_gene_list, Peri_DESeq2_gene_list)
Peri_venn <- ggVennDiagram(Peri_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Peri", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Peri_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

Endo_venn_list <- list(Endo_Limma_gene_list, Endo_EdgeR_gene_list, Endo_DESeq2_gene_list)
Endo_venn <- ggVennDiagram(Endo_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Endo", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Endo_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

unique_venn_list <- list(unique_Limma_genes, unique_EdgeR_genes, unique_DESeq2_genes)
unique_genes_venn <- ggVennDiagram(unique_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Unique Differentially Expressed Genes Identified for All Cell Types", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("unique_genes_M_MUT_and_WT_M_E18_WB_venn.pdf", device = "pdf", path = venn_dir)

# Show genes identified by all methods for cell types
Reduce(intersect, list(L2_3_IT_Limma_gene_list, L2_3_IT_DESeq2_gene_list, L2_3_IT_EdgeR_gene_list))
Reduce(intersect, list(L6_Limma_gene_list, L6_DESeq2_gene_list, L6_EdgeR_gene_list))
Reduce(intersect, list(Sst_Limma_gene_list, Sst_DESeq2_gene_list, Sst_EdgeR_gene_list))
Reduce(intersect, list(L5_Limma_gene_list, L5_DESeq2_gene_list, L5_EdgeR_gene_list))
Reduce(intersect, list(L4_Limma_gene_list, L4_DESeq2_gene_list, L4_EdgeR_gene_list))
Reduce(intersect, list(Pvalb_Limma_gene_list, Pvalb_DESeq2_gene_list, Pvalb_EdgeR_gene_list))
Reduce(intersect, list(Sncg_Limma_gene_list, Sncg_DESeq2_gene_list, Sncg_EdgeR_gene_list))
Reduce(intersect, list(Non_neuronal_Limma_gene_list, Non_neuronal_DESeq2_gene_list, Non_neuronal_EdgeR_gene_list))
Reduce(intersect, list(Oligo_Limma_gene_list, Oligo_DESeq2_gene_list, Oligo_EdgeR_gene_list))
Reduce(intersect, list(Vip_Limma_gene_list, Vip_DESeq2_gene_list, Vip_EdgeR_gene_list))
Reduce(intersect, list(Lamp5_Limma_gene_list, Lamp5_DESeq2_gene_list, Lamp5_EdgeR_gene_list))
Reduce(intersect, list(Astro_Limma_gene_list, Astro_DESeq2_gene_list, Astro_EdgeR_gene_list))
Reduce(intersect, list(Peri_Limma_gene_list, Peri_DESeq2_gene_list, Peri_EdgeR_gene_list))
Reduce(intersect, list(Endo_Limma_gene_list, Endo_DESeq2_gene_list, Endo_EdgeR_gene_list))

unique_genes <- Reduce(intersect, list(unique_Limma_genes, unique_DESeq2_genes, unique_EdgeR_genes))
for (unique_gene in unique_genes){
  print(unique_gene)
}
