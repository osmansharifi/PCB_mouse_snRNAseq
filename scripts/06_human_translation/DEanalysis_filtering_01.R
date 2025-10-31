#### Differential Expression Analysis of scRNA-seq Data - 01 - Create Design Matrix and Filtering Genes ####

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis", "scCustomize", "Seurat")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

load('/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/08_human_hdWGCNA/Rett_human_cortex.RData')
s.obj.name = "rett_P60_with_labels_proportions" #change this to the name of the Seurat object you started with

setwd(glue::glue("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis"))

s.obj <- human_rettcort
s.obj@meta.data$PCB_binary[grepl("4591_CTRL|662_CTRL", s.obj@meta.data$Samples)] <- NA
# Plot the graph, faceting by 'Group'
ggplot(s.obj@meta.data, aes(x = orig.ident, fill = PCB_binary)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("No" = "blue", "Yes" = "red")) +  # Custom colors for "No" and "Yes"
  labs(title = "PCB Binary Status by Sample", x = "Sample", y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  facet_wrap(~ Condition, scales = "free_x")+  # Create separate plots for CTRL and RTT
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = "bottom")+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 45, size = 18, face = 'bold', hjust = 0.4, vjust = .6, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 18, face = 'bold', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"), # Text size
    title = element_text(size = 18, face = "bold")) 
ggsave(glue::glue("PCB_presence_persample.pdf"), 
       width = 15,
       height = 12)

## Subset only the RTT samples to perfrom DEG analysis on
Idents(s.obj) <- "Samples"
RTT <- subset(x = s.obj, idents = c("1420_RTT", "1815_RTT", "3381_RTT", "4687_RTT"))
Idents(RTT) <- "predicted.class"
## creating design matrix 

expr_matrix = RTT@assays$RNA@counts

design = data.frame(cell_type = RTT$Cell_type,
                    sample_ID = RTT$orig.ident,
                    cell_cycle = RTT$cell.cycle,
                    Treatment = ifelse(grepl("Yes", RTT$orig.ident)=="TRUE", "Yes", "No"),
                    cell_ID = colnames(expr_matrix),
                    percent.mito = RTT$percent.mito)

design = design %>%
  dplyr::mutate_if(is.character, as.factor)

design$Treatment = factor(design$Treatment, levels=c("Yes", "No"))

## creating count matrices split by cell type

cell_types = as.factor(levels(design$cell_type))

expr_matrix_list = lapply(levels(cell_types), function(x) {
  expr_matrix[,design$cell_ID[which(design$cell_type==x)]]
})

names(expr_matrix_list) = as.character(cell_types)

expr_matrix_new = c(expr_matrix_list)

# create DGEList object

cell_types_all = names(expr_matrix_new)

DGEList = lapply(cell_types_all, function(x) {
  DGEList(expr_matrix_new[[x]])
})

names(DGEList) = cell_types_all

### Filtering DGEList by highly and lowly expressed genes. Highly expressed = at least 1 CPM in more than 25% of cells

DGEListCPM = lapply(cell_types_all, function(x){
  cpm(DGEList[[x]])
})

names(DGEListCPM) = cell_types_all

highly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})

names(highly_expr_genes) = cell_types_all

data = sapply(cell_types_all, function(x) {length(highly_expr_genes[[x]])})
data.frame1 = data.frame(cell_type = cell_types_all, num_highly_expressed_genes=data, cut_off = rep("25%", times=length(cell_types_all)))

highly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.50*ncol(DGEList[[x]]$counts))=="TRUE")
})

names(highly_expr_genes) = cell_types_all

data = sapply(cell_types_all, function(x) {length(highly_expr_genes[[x]])})
data.frame2 = data.frame(cell_type = cell_types_all, num_highly_expressed_genes=data, cut_off = rep("50%", times=length(cell_types_all)))

data.frame.comb = rbind(data.frame1, data.frame2)

pdf(file="Gene_Filtering_cutoffs_25_50.pdf", height=8.5, width=11)
ggplot(data=data.frame.comb, aes(x=cell_type, y=num_highly_expressed_genes)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  facet_grid(~cut_off) +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="NULL") +
  ggtitle("Genes with at least 1 CPM in more than either 25% or 50% of cells") 
dev.off()

# Using 25% cutoff
highly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})
names(highly_expr_genes) = cell_types_all

lowly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) < 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})
names(lowly_expr_genes) = cell_types_all

DGEList_high = lapply(cell_types_all, function(x){
  DGEList[[x]][highly_expr_genes[[x]],] %>%
    calcNormFactors()
})

DGEList_low = lapply(cell_types_all, function(x){
  DGEList[[x]][lowly_expr_genes[[x]],] %>%
    calcNormFactors()
})

names(DGEList_high) = cell_types_all
names(DGEList_low) = cell_types_all

save(DGEList, DGEList_high, design, cell_types_all, file="DEanalysis_01.RData")


