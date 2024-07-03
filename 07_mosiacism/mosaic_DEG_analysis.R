library(Seurat)
library(dplyr)
library(scCustomize)
library(patchwork)
library(glue)
library(limma)
library(ggplot2)
library(ggrepel)
library(edgeR)
##################
## Load samples ##
##################
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
load(glue('{base_path}/all.female.cortex.parsed.RData'))
mosaic.cortex <- subset(x = all.female.cortex, subset = Condition == 'MUTANT')
cluster <- subset(mosaic.cortex, idents = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5", "Oligo","Astro","Non-neuronal"))

###################################################################
## Limma DEG analysis of mutant vs wt cells from the het females ##
###################################################################
# Get expression info 
expr <- as.matrix(GetAssayData(cluster))
# Filter out genes that are 0 for every cell in this cluster
bad <- which(rowSums(expr) == 0)
expr <- expr[-bad,]
logcpm <- cpm(expr, prior.count=2, log=TRUE)
mm <- model.matrix(~0 + Mecp2_allele, data = cluster@meta.data)
y <- voom(expr, mm, plot = TRUE)
fit <- lmFit(y, mm)  
head(coef(fit)) # means in each sample for each gene
contr<- makeContrasts(c(Mecp2_alleleMecp2_MUT) - c(Mecp2_alleleMecp2_WT), levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrasts = contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
length(which(top.table$adj.P.Val < 0.05))
DEG<-as.matrix(topTable(tmp, sort.by = "adj.P.Val", n = 100))
summary(decideTests(tmp))

###########################
## Prepare and Visualize ##
###########################
top.table$Gene <- rownames(top.table)
top.table$diffexpressed <- 'NO'
top.table$diffexpressed[top.table$logFC > 0 & top.table$adj.P.Val < 0.05] <- 'UP'
top.table$diffexpressed[top.table$logFC < 0 & top.table$adj.P.Val < 0.05] <- 'DOWN'
top.table$diffexpressed[top.table$adj.P.Val > 0.05] <- 'Not Sig'
top.table$delabel <- NA
thresh = head(arrange(top.table, adj.P.Val), 10)$adj.P.Val[10]
top.table$delabel[top.table$adj.P.Val <=thresh] <-(top.table$Gene[top.table$adj.P.Val<=thresh])

# Volcano Plot
ggplot(data = top.table, aes(x = logFC, y = -log(adj.P.Val), col = diffexpressed, label = delabel))+
  geom_point(size=3)+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c('blue', 'black', 'red'))+
  theme(text = element_text(size=16)) +
theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
    labs(title = 'Mutant Cells VS WT Cells within the mosaic brains')
ggplot2::ggsave(glue("{base_path}/VolcanoPlot_MUTvsWT_withinfemales.tiff"),
                device = NULL,
                height = 8.5,
                width = 12)

#####################################
## Testing DEGs for each timepoint ##
#####################################
# Perform DEG analysis for each age group
age_groups <- unique(mosaic.cortex@meta.data$Age)
deg_results <- list()

for (age_group in age_groups) {
  cat("Performing DEG analysis for", age_group, "\n")
  
  # Subset cells based on age
  age_subset <- subset(mosaic.cortex, subset = Age == age_group)
  
  # Get expression info
  expr <- as.matrix(GetAssayData(age_subset))
  
  # Filter out genes that are 0 for every cell
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad, ]
  
  logcpm <- cpm(expr, prior.count = 2, log = TRUE)
  mm <- model.matrix(~0 + Mecp2_allele, data = age_subset@meta.data)
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Extract DEG results
  contrasts <- makeContrasts(c(Mecp2_alleleMecp2_MUT) - c(Mecp2_alleleMecp2_WT), levels = colnames(coef(fit)))
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
deg_results$P30$Timepoint <- 'P30'
deg_results$P60$Timepoint <- 'P60'
deg_results$P150$Timepoint <- 'P150'
deg_results$P30$DEG_Test <- 'MUT cells from HET vs WT cells from HET'
deg_results$P60$DEG_Test <- 'MUT cells from HET vs WT cells from HET'
deg_results$P150$DEG_Test <- 'MUT cells from HET vs WT cells from HET'
MUTvsWT_within_mosaic <- rbind(deg_results$P30, deg_results$P60, deg_results$P150)
write.csv(MUTvsWT_within_mosaic, file = glue('{base_path}/mutvswt_DEG.csv'))

top.table <- deg_results$P150
top.table$Gene <- rownames(top.table)
top.table$diffexpressed <- 'NO'
top.table$diffexpressed[top.table$logFC > 0 & top.table$adj.P.Val < 0.05] <- 'UP'
top.table$diffexpressed[top.table$logFC < 0 & top.table$adj.P.Val < 0.05] <- 'DOWN'
top.table$diffexpressed[top.table$adj.P.Val > 0.05] <- 'Not Sig'
top.table$delabel <- NA
thresh = head(arrange(top.table, adj.P.Val), 30)$adj.P.Val[30]
top.table$delabel[top.table$adj.P.Val <=thresh] <-(top.table$Gene[top.table$adj.P.Val<=thresh])

# Volcano Plot
ggplot(data = top.table, aes(x = logFC, y = -log(adj.P.Val), col = diffexpressed, label = delabel))+
  geom_point(size=3)+
  theme_minimal()+
  geom_text_repel(min.segment.length = 0)+
  scale_color_manual(values = c('blue', 'black', 'red'))+
  theme(text = element_text(size=16)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  labs(title = 'Mutant Cells VS WT Cells within the P150 mosaic brains')
ggplot2::ggsave(glue("{base_path}/VolcanoPlot_MUTvsWT_withinP150females.tiff"),
                device = NULL,
                height = 8.5,
                width = 12)

###############################
## Testing DEGs for WT vs WT ##
###############################
# Perform DEG analysis between the WT cells from the WT mouse and WT cells from the mosaic brains
cell_nonautonomous <- subset(x = all.female.cortex, subset = Mecp2_allele == 'Mecp2_WT')
WT_from_MUT = Cells(cell_nonautonomous)[which(cell_nonautonomous$Condition == "MUTANT")]
WT_from_WT = Cells(cell_nonautonomous)[which(cell_nonautonomous$Condition == "WT")]
slct_WT_from_MUT = sample(WT_from_MUT, size = 607)
slct_WT_from_WT = sample(WT_from_WT, size = 607)
subset_cell_nonautonomous = subset(cell_nonautonomous, cells = c(slct_WT_from_MUT, slct_WT_from_WT))

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

top.table <- deg_results$P60
top.table$Gene <- rownames(top.table)
top.table$diffexpressed <- 'NO'
top.table$diffexpressed[top.table$logFC > 0 & top.table$adj.P.Val < 0.05] <- 'UP'
top.table$diffexpressed[top.table$logFC < 0 & top.table$adj.P.Val < 0.05] <- 'DOWN'
top.table$diffexpressed[top.table$adj.P.Val > 0.05] <- 'Not Sig'
top.table$delabel <- NA
thresh = head(arrange(top.table, adj.P.Val), 10)$adj.P.Val[10]
top.table$delabel[top.table$adj.P.Val <=thresh] <-(top.table$Gene[top.table$adj.P.Val<=thresh])

# Volcano Plot
ggplot(data = top.table, aes(x = logFC, y = -log(adj.P.Val), col = diffexpressed, label = delabel))+
  geom_point(size=3)+
  theme_minimal()+
  geom_text_repel(max.overlaps = Inf)+
  scale_color_manual(values = c('blue', 'black', 'red'))+
  theme(text = element_text(size=16)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  labs(title = 'WT cells from WT P60 Female vs WT cells from P60 mosaic cortex')
 ggplot2::ggsave(glue("{base_path}/VolcanoPlot_WTvsWT_P60females.tiff"),
                device = NULL,
                height = 8.5,
                width = 12)

 #Write csv files
WTvsWT <- rbind(deg_results$P30, deg_results$P60, deg_results$P150)
 ###########################################
 ## Venn diagram of the overlapping genes ##
 ###########################################
 library(VennDiagram)

 # Extract significant DEGs
 sig_genes_P30 <- rownames(deg_results$P30)[deg_results$P30$adj.P.Val < 0.05]
 sig_genes_P60 <- rownames(deg_results$P60)[deg_results$P60$adj.P.Val < 0.05]
 sig_genes_P150 <- rownames(deg_results$P150)[deg_results$P150$adj.P.Val < 0.05]
 intersection_all2 <- intersect(sig_genes_P150,sig_genes_P30)
 intersection_all3 <- intersect(intersection_all2,sig_genes_P60)
 # Create a Venn diagram
 venn.diagram(
   x = list(
     P30 = sig_genes_P30,
     P60 = sig_genes_P60,
     P150 = sig_genes_P150
   ),
   category.names = c("P30", "P60", "P150"),
   main = 'DEGs from WT cells from WT females and WT cells from mosaic females ',
   filename = glue("{base_path}/venn_test.tiff"),
   col = c('red', 'green', 'blue'),
   fill = c('red', 'green', 'blue'),
   cat.cex = 1.2,
   cat.fontface = "bold",
   euler.d = TRUE,
   disable.logging = TRUE,
   hyper.test = TRUE
 )

DEGs = top.table
# Top DEGs

DEGs <- fit %>%
  topTable(sort.by = "P", n = Inf) %>%
  tibble::rownames_to_column() %>%
  tibble::as_tibble() %>%
  dplyr::rename(SYMBOL = rowname) %>%
  dplyr::mutate(FC = dplyr::case_when(logFC >0 ~ 2^logFC,
                                      logFC <0 ~ -1/(2^logFC))) %>%
  dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %T>%
  openxlsx::write.xlsx(file=glue::glue("{i}/DEGs.xlsx")) %>%
  dplyr::filter(P.Value < 0.05) %T>%
  openxlsx::write.xlsx(file=glue::glue("{i}/sig_DEGs.xlsx"))

print(glue::glue("GO and Pathway analysis of {i} cells"))

enrichR:::.onAttach()
source(glue('{base_path}/GO_ploting_functions.R'))
tryCatch({
  DEGs %>% 
    dplyr::select(Gene) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Molecular_Function_2018",
                       "GO_Cellular_Component_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
    openxlsx::write.xlsx(file=glue::glue("{base_path}/enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{base_path}/rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("{base_path}/enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
  })
print(glue::glue("The pipeline has finished for samples"))

