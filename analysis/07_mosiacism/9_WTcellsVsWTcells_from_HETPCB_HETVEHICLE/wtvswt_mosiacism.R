########################################################
## Create broad categories and run mosiacism analysis ##
########################################################
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
##################
## Load samples ##
##################
base_path <- '/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/'
load(glue('{base_path}/PEBBLES_parsed.RData'))
mosaic.cortex <- subset(x = PEBBLES_soupx, subset = Genotype == 'HET')

# Perform DEG analysis between the WT cells from the WT mouse and WT cells from the mosaic brains
cell_nonautonomous <- subset(x = mosaic.cortex, subset = Mecp2_allele == 'WT_Mecp2')
WT_from_PCB = Cells(cell_nonautonomous)[which(cell_nonautonomous$Treatment == "PCB")]
WT_from_VEHICLE = Cells(cell_nonautonomous)[which(cell_nonautonomous$Treatment == "VEHICLE")]
slct_WT_from_PCB = sample(WT_from_PCB, size = 199)
slct_WT_from_VEHICLE = sample(WT_from_VEHICLE, size = 199)
subset_cell_nonautonomous = subset(cell_nonautonomous, cells = c(slct_WT_from_PCB, slct_WT_from_VEHICLE))
celltypes <- unique(subset_cell_nonautonomous@meta.data$broad_class)
deg_results <- list()

for (celltype in celltypes) {
  cat("Performing DEG analysis for", celltype, "\n")
  
  # Subset cells based on broad_class
  broad_class_subset <- subset(subset_cell_nonautonomous, subset = broad_class == celltype)
  
  # Get expression info
  expr <- as.matrix(GetAssayData(broad_class_subset))
  
  # Filter out genes that are 0 for every cell
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad, ]
  
  logcpm <- cpm(expr, prior.count = 2, log = TRUE)
  mm <- model.matrix(~0 + Treatment, data = broad_class_subset@meta.data)
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Extract DEG results
  contrasts <- makeContrasts(c(TreatmentPCB) - c(TreatmentVEHICLE), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contrasts)
  tmp <- eBayes(tmp)
  top_table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
  
  # Store DEG results for this broad_class group
  deg_results[[celltype]] <- top_table
}

# Access DEG results for each broad_class group
for (celltype in celltypes) {
  cat("DEG analysis results for", celltype, "\n")
  
  deg_table <- deg_results[[celltype]]
  num_degs <- length(which(deg_table$adj.P.Val < 0.05))
  
  print(num_degs)
  
  # Additional analysis if needed
  # summary(decideTests(tmp))
}

top.table <- deg_results$Glutamatergic
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
directory_path = glue('{base_path}9_WTcellsVsWTcells_from_HETPCB_HETVEHICLE')
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
ggplot2::ggsave(glue("{directory_path}/Glutamatergic_Volcano_WTvsWT_HetPCBvshHetVEH.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

#Write csv files
deg_results$Glutamatergic$Cell_type <- "Glutamatergic"
deg_results$GABAergic$Cell_type <- "GABAergic"
deg_results$`Non-neuronal`$Cell_type <- "Non-neuronal"
WTvsWT <- rbind(deg_results$Glutamatergic, deg_results$GABAergic, deg_results$`Non-neuronal`)
WTvsWT$SYMBOL <- rownames(WTvsWT)
WTvsWT_significant <- subset(WTvsWT, adj.P.Val <= 0.05)
write.csv(WTvsWT_significant, file = glue('{directory_path}/sig_WTcellsVsWTcells_from_HETPCB_HETVEHICLE_DEGs.csv'))
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
pdf(glue("{directory_path}/venn_sig_WTcellsVsWTcells_from_HETPCB_HETVEHICLE_DEGs.pdf"))
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

######################
## Pathway analysis ##
######################
DEGs = top.table
# Top DEGs

DEGs <- deg_results$Glutamatergic %>%
  tibble::rownames_to_column() %>%
  tibble::as_tibble() %>%
  dplyr::rename(SYMBOL = rowname) %>%
  dplyr::mutate(FC = dplyr::case_when(logFC >0 ~ 2^logFC,
                                      logFC <0 ~ -1/(2^logFC))) %>%
  dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %T>%
  openxlsx::write.xlsx(file=glue::glue("{directory_path}/DEGs.xlsx")) %>%
  dplyr::filter(P.Value < 0.05) %T>%
  openxlsx::write.xlsx(file=glue::glue("{directory_path}/sig_DEGs.xlsx"))

print(glue::glue("GO and Pathway analysis of {i} cells"))

enrichR:::.onAttach()
source('/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/scripts/09_mosiacism_analysis/GO_ploting_functions.R')
tryCatch({
  DEGs %>% 
    dplyr::select(SYMBOL) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Molecular_Function_2023",
                       "GO_Cellular_Component_2023",
                       "KEGG_2019_Mouse",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
    openxlsx::write.xlsx(file=glue::glue("{directory_path}/Glutamatergic_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Mm.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{directory_path}/Glutamatergic_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("{directory_path}/GABAergic_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
  })
print(glue::glue("The pipeline has finished for samples"))
