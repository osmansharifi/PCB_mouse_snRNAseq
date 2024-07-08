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
##################
## Load samples ##
##################
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
load(glue('{base_path}/all.female.cortex.parsed.RData'))
mosaic.cortex <- subset(x = all.female.cortex, subset = Condition == 'MUTANT')
cluster <- subset(mosaic.cortex, idents = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5", "Oligo","Astro","Non-neuronal"))

# Create a new column called broad_class based on celltype.call
all.female.cortex$broad_class <- ifelse(
  all.female.cortex$celltype.call %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Vip"), 
  "GABAergic", 
  ifelse(
    all.female.cortex$celltype.call %in% c("L2_3_IT", "L4", "L5", "L6"), 
    "Glutamatergic", 
    ifelse(
      all.female.cortex$celltype.call %in% c("Astro", "Non-neuronal", "Oligo"), 
      "Non-neuronal", 
      "Other"
    )
  )
)

# Perform DEG analysis between the WT cells from the WT mouse and WT cells from the mosaic brains
cell_nonautonomous <- subset(x = all.female.cortex, subset = Mecp2_allele == 'Mecp2_WT')
WT_from_MUT = Cells(cell_nonautonomous)[which(cell_nonautonomous$Condition == "MUTANT")]
WT_from_WT = Cells(cell_nonautonomous)[which(cell_nonautonomous$Condition == "WT")]
slct_WT_from_MUT = sample(WT_from_MUT, size = 607)
slct_WT_from_WT = sample(WT_from_WT, size = 607)
subset_cell_nonautonomous = subset(cell_nonautonomous, cells = c(slct_WT_from_MUT, slct_WT_from_WT))
subset_cell_nonautonomous <- subset(x = subset_cell_nonautonomous, subset = broad_class == 'Glutamatergic')
age_groups <- unique(subset_cell_nonautonomous@meta.data$Age)
deg_results <- list()

for (age_group in age_groups) {
  cat("Performing DEG analysis for", age_group, "\n")
  
  # Subset cells based on age
  age_subset <- subset(subset_cell_nonautonomous, subset = Age == age_group)
  
  # Get expression info
  expr <- as.matrix(GetAssayData(age_subset))
  
  # Filter out genes that are 0 for every cell
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad, ]
  
  logcpm <- cpm(expr, prior.count = 2, log = TRUE)
  mm <- model.matrix(~0 + Condition, data = age_subset@meta.data)
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Extract DEG results
  contrasts <- makeContrasts(c(ConditionMUTANT) - c(ConditionWT), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contrasts)
  tmp <- eBayes(tmp)
  top_table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
  
  # Store DEG results for this age group
  deg_results[[age_group]] <- top_table
}

# Access DEG results for each age group
for (age_group in age_groups) {
  cat("DEG analysis results for", age_group, "\n")
  
  deg_table <- deg_results[[age_group]]
  num_degs <- length(which(deg_table$adj.P.Val < 0.05))
  
  print(num_degs)
  
  # Additional analysis if needed
  # summary(decideTests(tmp))
}

top.table <- deg_results$P30
top.table$Gene <- rownames(top.table)
# Add necessary columns to the data frame
DEGs$Gene <- DEGs$SYMBOL
DEGs$diffexpressed <- 'NO'
DEGs$diffexpressed[DEGs$logFC > 0 & DEGs$adj.P.Val < 0.05] <- 'UP'
DEGs$diffexpressed[DEGs$logFC < 0 & DEGs$adj.P.Val < 0.05] <- 'DOWN'
DEGs$diffexpressed[DEGs$adj.P.Val > 0.05] <- 'Not Sig'

# Get the number of significant upregulated and downregulated genes
num_upregulated <- sum(DEGs$logFC > 0 & DEGs$adj.P.Val < 0.05)
num_downregulated <- sum(DEGs$logFC < 0 & DEGs$adj.P.Val < 0.05)

# Get the top 10 upregulated genes
top_upregulated_genes <- DEGs %>% arrange(desc(logFC)) %>% head(5)

# Get the top 10 downregulated genes
top_downregulated_genes <- DEGs %>% arrange(logFC) %>% head(5)

# Create a column for label based on top genes
DEGs$delabel <- NA
DEGs$delabel[DEGs$Gene %in% top_upregulated_genes$Gene] <- top_upregulated_genes$Gene
DEGs$delabel[DEGs$Gene %in% top_downregulated_genes$Gene] <- top_downregulated_genes$Gene

# Get the directory name from the directory path
dir_name <- basename(directory_path)

# Volcano Plot
gg <- ggplot(data = DEGs, aes(x = logFC, y = -log(adj.P.Val), col = diffexpressed, label = delabel)) +
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
       subtitle = paste("Upregulated:", num_upregulated, " | Downregulated:", num_downregulated)) +  # Add subtitle with counts
  
  # Set x-axis limits from -2 to 2
  xlim(c(-0.5, 0.75))+
  ylim(c(0, 70))
# Save the plot in the current directory with the name "<dir_name>_volcano.pdf"
ggsave(plot = gg,
       filename = file.path(directory_path, paste0(dir_name, "_volcano.pdf")),
       height = 8.5,
       width = 12,
       device = NULL)
  labs(title = 'Glutamatergic WT cells from WT P150 females vs Glutamatergic WT cells from P150 mosaic')
ggplot2::ggsave(glue("{base_path}/broad_group_analysis/Vol_glut_WTvsWT_P150females.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

ggplot2::ggsave(glue("{base_path}/broad_group_analysis/celltype_UMAP.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)
#Write csv files
deg_results$P30$Time_point <- "P30"
deg_results$P60$Time_point <- "P60"
deg_results$P150$Time_point <- "P150"
WTvsWT <- rbind(deg_results$P30, deg_results$P60, deg_results$P150)
WTvsWT$SYMBOL <- rownames(WTvsWT)
write.csv(WTvsWT, file = glue('{base_path}WTvsWT_DEGs.csv'))
###########################################
## Venn diagram of the overlapping genes ##
###########################################

# Extract significant DEGs
sig_genes_P30 <- rownames(deg_results$P30)[deg_results$P30$adj.P.Val < 0.05]
sig_genes_P60 <- rownames(deg_results$P60)[deg_results$P60$adj.P.Val < 0.05]
sig_genes_P150 <- rownames(deg_results$P150)[deg_results$P150$adj.P.Val < 0.05]
intersection_all2 <- intersect(sig_genes_P150,sig_genes_P30)
intersection_all3 <- intersect(intersection_all2,sig_genes_P60)
intersection_all4 <- intersect(sig_genes_P30,sig_genes_P60)
# Create a Venn diagram
pdf(glue("{base_path}/broad_group_analysis/venn_Glutamatergic.pdf"))
temp <- venn.diagram(
  x = list(
    P30 = sig_genes_P30,
    P60 = sig_genes_P60,
    P150 = sig_genes_P150
  ),
  category.names = c("P30", "P60", "P150"),
  main = 'Glutamatergic DEGs from WT cells from WT females and WT cells from mosaic females ',
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
max_length <- max(length(sig_genes_P30), length(sig_genes_P60), length(sig_genes_P150))
sig_genes_P30 <- c(sig_genes_P30, rep(NA, max_length - length(sig_genes_P30)))
sig_genes_P60 <- c(sig_genes_P60, rep(NA, max_length - length(sig_genes_P60)))
sig_genes_P150 <- c(sig_genes_P150, rep(NA, max_length - length(sig_genes_P150)))

# Combine into a dataframe
combined_df <- data.frame(sig_genes_P30, sig_genes_P60, sig_genes_P150)

go.obj <- newGeneOverlap(combined_df$sig_genes_P30,
                         combined_df$sig_genes_P60,
                         combined_df$sig_genes_P150,
                         genome.size = 22000)
go.obj <- testGeneOverlap(go.obj)
getPval(go.obj)
getOddsRatio(go.obj)
getJaccard(go.obj)
getContbl(go.obj)
print(go.obj)
write.csv(as.data.frame(go.obj), file = glue('{base_path}/broad_group_analysis/geneoverlap_gaba_3timepoints.txt'))
DEGs = top.table
# Top DEGs

DEGs <- deg_results$P60 %>%
  tibble::rownames_to_column() %>%
  tibble::as_tibble() %>%
  dplyr::rename(SYMBOL = rowname) %>%
  dplyr::mutate(FC = dplyr::case_when(logFC >0 ~ 2^logFC,
                                      logFC <0 ~ -1/(2^logFC))) %>%
  dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %T>%
  openxlsx::write.xlsx(file=glue::glue("{base_path}/broad_group_analysis/DEGs.xlsx")) %>%
  dplyr::filter(P.Value < 0.05) %T>%
  openxlsx::write.xlsx(file=glue::glue("{base_path}/broad_group_analysis/sig_DEGs.xlsx"))

print(glue::glue("GO and Pathway analysis of {i} cells"))

enrichR:::.onAttach()
source(glue('{base_path}/GO_ploting_functions.R'))
tryCatch({
  DEGs %>% 
    dplyr::select(SYMBOL) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Molecular_Function_2018",
                       "GO_Cellular_Component_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
    openxlsx::write.xlsx(file=glue::glue("{base_path}/broad_group_analysis/GABA_P60_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Mm.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{base_path}/broad_group_analysis/GABA_P60_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("{base_path}/broad_group_analysis/GABA_P60_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
  })
print(glue::glue("The pipeline has finished for samples"))
