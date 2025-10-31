#### Differential Expression Analysis of scRNA-seq Data - 01 - Create Design Matrix and Filtering Genes ####

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

setwd(glue::glue("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/"))
load("PEBBLES_soupx_labeled.RData")
s.obj <- subset(PEBBLES_soupx, subset = Treatment == "VEHICLE")

## creating design matrix 

expr_matrix = s.obj@assays$RNA@counts

design = data.frame(cell_type = s.obj$predicted.id,
                    #activation.status = ifelse(grepl("activated", s.obj$cell.type)=="TRUE", "activated", "unactivated"),
                    sample_ID = s.obj$Samples,
                    #cell_cycle = s.obj$cell.cycle,
                    cell_ID = colnames(expr_matrix),
                    genotype = s.obj$Genotype,
                    percent.mito = s.obj$percent.mito,
                    treatment = s.obj$Treatment)

design = design %>%
  dplyr::mutate_if(is.character, as.factor)

design$treatment = factor(design$treatment)
design$sample_ID = factor(design$sample_ID)
design$genotype = factor(design$genotype, levels = c("HET", "WT"))
## creating count matrices split by cell type

cell_types = as.factor(levels(design$cell_type))

expr_matrix_list = lapply(levels(cell_types), function(x) {
  expr_matrix[,design$cell_ID[which(design$cell_type==x)]]
})

names(expr_matrix_list) = as.character(cell_types)

# create DGEList object

cell_types_all = names(expr_matrix_list)

DGEList = lapply(cell_types_all, function(x) {
  DGEList(expr_matrix_list[[x]])
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

#lowly_expr_genes = lapply(cell_types_all, function(x){
  #which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) < 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
#})
#names(lowly_expr_genes) = cell_types_all

DGEList_high = lapply(cell_types_all, function(x){
  DGEList[[x]][highly_expr_genes[[x]],] %>%
    calcNormFactors()
})

#DGEList_low = lapply(cell_types_all, function(x){
#  DGEList[[x]][lowly_expr_genes[[x]],] %>%
#    calcNormFactors()
#})

names(DGEList_high) = cell_types_all
#names(DGEList_low) = cell_types_all

save(DGEList, DGEList_high, design, cell_types_all, file="DEanalysis_WTvsHET_in_VEHICLE_01.RData")


