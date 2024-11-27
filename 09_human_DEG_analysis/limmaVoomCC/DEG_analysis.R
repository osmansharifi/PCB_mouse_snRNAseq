####################
## Load libraries ##
####################
library(magrittr)
library(VennDiagram)
library(grDevices)
library(dplyr)
library(Seurat)
library(glue)
library(scCustomize)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(GeneOverlap)
library(enrichR)
library(RColorBrewer)
##################
## Load samples ##
##################
base_path <- '/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/'
load('/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/08_human_hdWGCNA/Rett_human_cortex.RData')

s.obj <- subset(x = human_rettcort, subset = Condition == 'RTT')
s.obj <- subset(s.obj, subset = PCB_binary != 'NA')

# Perform DEG analysis between the WT cells from the WT mouse and WT cells from the mosaic brains
#cell_nonautonomous <- subset(x = s.obj, subset = Mecp2_allele == 'WT_Mecp2')
#WT_from_PCB = Cells(cell_nonautonomous)[which(cell_nonautonomous$Treatment == "PCB")]
#WT_from_VEHICLE = Cells(cell_nonautonomous)[which(cell_nonautonomous$Treatment == "VEHICLE")]
#slct_WT_from_PCB = sample(WT_from_PCB, size = 199)
#slct_WT_from_VEHICLE = sample(WT_from_VEHICLE, size = 199)
#subset_cell_nonautonomous = subset(cell_nonautonomous, cells = c(slct_WT_from_PCB, slct_WT_from_VEHICLE))
celltypes <- unique(s.obj@meta.data$predicted.class)
deg_results <- list()  # To store results for each cell type

for (celltype in celltypes) {
  cat("Performing DEG analysis for", celltype, "\n")
  
  # Subset cells based on Cell_type
  cell_type_subset <- subset(s.obj, subset = predicted.class == celltype)
  
  # Get expression data
  expr <- as.matrix(GetAssayData(cell_type_subset))
  
  # Filter out genes with zero counts in all cells
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad, ]
  
  # Check dimensions
  cat("Number of cells in expression matrix:", ncol(expr), "\n")
  cat("Number of cells in metadata:", nrow(cell_type_subset@meta.data), "\n")
  
  if (!all(colnames(expr) == rownames(cell_type_subset@meta.data))) {
    stop("Mismatch between expression matrix and metadata cell names.")
  }
  
  # Compute log-transformed CPM
  logcpm <- edgeR::cpm(expr, prior.count = 2, log = TRUE)
  
  # Create the design matrix
  mm <- model.matrix(~0 + PCB_binary, data = cell_type_subset@meta.data)
  
  # Ensure design matrix matches expression data dimensions
  cat("Design matrix dimensions:", dim(mm), "\n")
  if (nrow(mm) != ncol(expr)) {
    stop("Mismatch between design matrix and expression data.")
  }
  
  # Run voom normalization and linear model fitting
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Perform contrast fitting and eBayes
  contrasts <- makeContrasts(PCB_binaryYes - PCB_binaryNo, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contrasts)
  tmp <- eBayes(tmp)
  
  # Extract top DEG results
  top_table <- topTable(tmp, sort.by = "P", n = Inf)
  
  # Add Cell_type information
  top_table$Cell_type <- celltype
  
  # Store the DEG results for this cell type in the list
  deg_results[[celltype]] <- top_table
}

# Access DEG results for each predicted.class group
for (celltype in celltypes) {
  cat("DEG analysis results for", celltype, "\n")
  
  deg_table <- deg_results[[celltype]]
  num_degs <- length(which(deg_table$adj.P.Val < 0.05))
  
  print(num_degs)
}

top.table <- deg_results$GABAergic
# Subset top_table for rows where Cell_type is "Non-Neuronal"
#top.table <- top_table[top_table$Cell_type == `Non-Neuronal`, ]
top.table$Gene <- rownames(top.table)
# Add necessary columns to the data frame
top.table$diffexpressed <- 'NO'
top.table$diffexpressed[top.table$logFC > 0 & top.table$adj.P.Val < 0.05] <- 'UP'
top.table$diffexpressed[top.table$logFC < 0 & top.table$adj.P.Val < 0.05] <- 'DOWN'
top.table$diffexpressed[top.table$adj.P.Val > 0.05] <- 'Not Sig'

# Get the number of significant upregulated and downregulated genes
num_upregulated <- sum(top.table$logFC > 0 & top.table$adj.P.Val < 0.05)
num_downregulated <- sum(top.table$logFC < 0 & top.table$adj.P.Val < 0.05)

# Get the top 10 upregulated genes
top_upregulated_genes <- top.table %>% arrange(desc(logFC)) %>% head(5)

# Get the top 10 downregulated genes
top_downregulated_genes <- top.table %>% arrange(logFC) %>% head(5)

# Create a column for label based on top genes
top.table$delabel <- NA
top.table$delabel[top.table$Gene %in% top_upregulated_genes$Gene] <- top_upregulated_genes$Gene
top.table$delabel[top.table$Gene %in% top_downregulated_genes$Gene] <- top_downregulated_genes$Gene

# Get the directory name from the directory path
directory_path = glue('{base_path}/LimmaVoomCC')
dir_name <- basename(directory_path)

# Volcano Plot
ggplot(data = top.table, aes(x = logFC, y = -log(adj.P.Val), col = diffexpressed, label = delabel)) +
  geom_point(size=2) +
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values = c('blue', 'black', 'red')) +
  theme(
    text = element_text(size=16),
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 0, size = 16, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 16, face = 'bold'),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")
  ) +
  labs(title = paste("Volcano Plot -", dir_name),  # Update the plot title
       subtitle = paste("Upregulated:", num_upregulated, " | Downregulated:", num_downregulated))   # Add subtitle with counts
ggplot2::ggsave(glue("{directory_path}/Nonneuronal_Volcano_{dir_name}.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

#Write csv files
deg_results$Glutamatergic$Cell_type <- "Glutamatergic"
deg_results$GABAergic$Cell_type <- "GABAergic"
deg_results$`Non-neuronal`$Cell_type <- "Non_Neuronal"
all_degs <- rbind(deg_results$Glutamatergic, deg_results$GABAergic, deg_results$`Non-Neuronal`)
all_degs$SYMBOL <- rownames(all_degs)
sig_all_degs <- all_degs[all_degs$adj.P.Val <= 0.05, ]
Human_rett_PCBvs_noPCB_sig_DEGs = sig_all_degs
write.csv(Human_rett_PCBvs_noPCB_sig_DEGs, file = glue('{directory_path}/Human_rett_PCBvs_noPCB_sig_DEGs.csv'))
###########################################
## Venn diagram of the overlapping genes ##
###########################################

# Extract significant DEGs
sig_genes_Glutamatergic <- rownames(deg_results$Glutamatergic)[deg_results$Glutamatergic$adj.P.Val < 0.05]
sig_genes_GABAergic <- rownames(deg_results$GABAergic)[deg_results$GABAergic$adj.P.Val < 0.05]
sig_genes_Non_neuronal <- rownames(deg_results$`Non-Neuronal`)[deg_results$`Non-Neuronal`$adj.P.Val < 0.05]
intersection_all2 <- intersect(sig_genes_Non_neuronal,sig_genes_Glutamatergic)
intersection_all3 <- intersect(intersection_all2,sig_genes_GABAergic)
intersection_all4 <- intersect(sig_genes_Glutamatergic,sig_genes_GABAergic)
# Create a Venn diagram
pdf(glue("{directory_path}/venn_sig_{dir_name}.pdf"))
temp <- venn.diagram(
  x = list(
    Glutamatergic = sig_genes_Glutamatergic,
    GABAergic = sig_genes_GABAergic,
    Non_neuronal = sig_genes_Non_neuronal
  ),
  category.names = c("Glutamatergic", "GABAergic", "Non-neuronal"),
  main = 'sig DEGs human Rett cortex PCB vs no PCB',
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

######################
## Pathway analysis ##
######################
DEGs = deg_results$GABAergic
# Ensure SYMBOL is a column (convert row names to SYMBOL)
DEGs <- DEGs %>%
  tibble::rownames_to_column("SYMBOL") %>%
  tibble::as_tibble() %>%
  dplyr::mutate(FC = dplyr::case_when(
    logFC > 0 ~ 2^logFC,
    logFC < 0 ~ -1/(2^logFC)
  )) %>%
  dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %>%
  arrange(desc(abs(logFC)))

DEGs <- DEGs %>%
  head(500)
# Sort by P.Value in ascending order
DEGs <- DEGs %>%
  arrange(P.Value)

# Write all DEGs to an Excel file
openxlsx::write.xlsx(DEGs, file = glue("{directory_path}/DEGs.xlsx"))

tryCatch({
  enriched_results <- enrichR::enrichr(DEGs$SYMBOL, c("GO_Biological_Process_2021",
                                         "GO_Molecular_Function_2021",
                                         "GO_Cellular_Component_2021",
                                         "KEGG_2021_Human",
                                         "Panther_2016",
                                         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"))
  
  print(str(enriched_results))  # Check the structure of the enrichment results
  enriched_results %>%
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
    openxlsx::write.xlsx(file=glue::glue("{directory_path}/GABAergic_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Hs.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{directory_path}/GABAergic_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("{directory_path}/GABAergic_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10)
}, error = function(error_condition) {
  print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
  print(error_condition) # Print the error for debugging
})

####
#Plot Pathways
####
library(openxlsx)
#load files
glut <- read.xlsx("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/Glutamatergic_enrichr.xlsx", sheet = 4)
gaba <- read.xlsx("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/GABAergic_enrichr.xlsx", sheet = 4)
NN <- read.xlsx("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/Non_neuronal_enrichr.xlsx", sheet = 4)

# Extract top 10 rows from each dataframe
df1_top <- head(glut, 10)
df2_top <- head(gaba, 10)
df3_top <- head(NN, 10)

# Combine the data frames into one
combined_df <- bind_rows(
  df1_top %>% mutate(Cell_Type = "Glutamatergic"),
  df2_top %>% mutate(Cell_Type = "GABAergic"),
  df3_top %>% mutate(Cell_Type = "Non-Neuronal")
)

# Plot using ggplot
ggplot(combined_df, aes(x = Term, y = Cell_Type, color = P.value)) +
  geom_point(size = 8) +
  scale_color_gradientn(name = "P.value", colors = brewer.pal(n = 11, name = "RdBu")) +
  scale_size_continuous(range = c(3, 8)) +
  scale_y_discrete(name = "Cell Type") +
  scale_x_discrete(name = "Terms") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'KEGG Terms') +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 90, size = 18, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 18, face = 'bold', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    # Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"), # Text size
    title = element_text(size = 18, face = "bold")
  ) +
  coord_flip()
ggplot2::ggsave(glue("{directory_path}/top10KEGG_human.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)
