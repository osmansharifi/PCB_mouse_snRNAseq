
  title: "Female_PCB_Cortex"
author: "Osman Sharifi"

library(Seurat)
library(biomaRt)
library(scCustomize)
library(ggplot2)

### The percentage of reads that map to the mitochondrial genome
# We use the set of all genes, in mouse these genes can be identified as those that begin with 'mt'.
#load data
load("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_soupx.RData")
PEBBLES_soupx$percent.mito <- PercentageFeatureSet(PEBBLES_soupx, pattern = "^mt-")
mito.genes <- grep("^mt-", rownames(PEBBLES_soupx), value = T)
percent.mito <- Matrix::colSums(GetAssayData(PEBBLES_soupx, slot = "counts")[mito.genes, ]) / Matrix::colSums(GetAssayData(PEBBLES_soupx, slot = "counts"))

## Calculate cell cycle using scran, add to meta data
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
# Convert to matrix for use in cycle
mat <- as.matrix(GetAssayData(PEBBLES_soupx))

# Convert rownames to ENSEMBL IDs, Using biomaRt
ensembl<- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
anno <- getBM( values=rownames(mat), attributes=c("mgi_symbol","ensembl_gene_id") , filters= "mgi_symbol"  ,mart=ensembl)
ord <- match(rownames(mat), anno$mgi_symbol) # use anno$mgi_symbol if via biomaRt
rownames(mat) <- anno$ensembl_gene_id[ord] # use anno$ensembl_gene_id if via biomaRt
drop <- which(is.na(rownames(mat)))
mat <- mat[-drop,]
cycles <- scran::cyclone(mat, pairs=mm.pairs)
tmp <- data.frame(cell.cycle = cycles$phases)
rownames(tmp) <- colnames(mat)
PEBBLES_soupx <- AddMetaData(PEBBLES_soupx, tmp)
Idents(PEBBLES_soupx) <- 'Samples'
do.call("cbind", tapply(PEBBLES_soupx$nFeature_RNA, Idents(PEBBLES_soupx),quantile,probs=seq(0,1,0.05)))
RidgePlot(PEBBLES_soupx, features="nFeature_RNA")
RidgePlot(PEBBLES_soupx, features="nCount_RNA")
VlnPlot(PEBBLES_soupx, features="percent.mito")
plot(sort(Matrix::rowSums(GetAssayData(PEBBLES_soupx) >= 3)) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")

#Normalize data
PEBBLES_soupx <- NormalizeData(
  object = PEBBLES_soupx,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

#Scale data
PEBBLES_soupx <- ScaleData(object = PEBBLES_soupx,verbose = TRUE, vars.to.regress = c("cell.cycle", "percent.mito"))

#Unfiltered plots
QC_Plot_UMIvsGene(seurat_object = PEBBLES_soupx, low_cutoff_gene = 600, high_cutoff_gene = 5000, low_cutoff_UMI = 500,
                  high_cutoff_UMI = 20000) +labs(x = "UMIs per nucleus", y = "Genes per nucleus" )
ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/prefilter_UMIvsGENE.pdf", 
       device = NULL,
       height = 8.5,
       width = 12)
QC_Plot_GenevsFeature(seurat_object = PEBBLES_soupx, feature1 = "percent.mito", low_cutoff_gene = 600,
                      high_cutoff_gene = 5000, high_cutoff_feature = 5)+labs(x = "Percent.mito per nucleus", y = "Genes per nucleus")
ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/prefilter_genevsmito.pdf", 
       device = NULL,
       height = 8.5,
       width = 12)

# Filter data 
PEBBLES_soupx <- subset(x = PEBBLES_soupx, 
                          subset= (nCount_RNA >= 500) & (nCount_RNA <= 20000) & 
                            (nFeature_RNA >= 600) & (nFeature_RNA <= 5000) &
                            (percent.mito < 0.5))
# Find variable genes
PEBBLES_soupx <- FindVariableFeatures(
  object = PEBBLES_soupx,
  selection.method = "vst", verbose = TRUE)
VariableFeaturePlot_scCustom(seurat_object = PEBBLES_soupx, num_features = 20, repel = TRUE,
                             y_axis_log = TRUE)
ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/variable_genes.pdf", 
       device = NULL,
       height = 8.5,
       width = 12)

#Run PCA
PEBBLES_soupx <- RunPCA(object = PEBBLES_soupx, seed.use = 1234)
DimPlot_scCustom(seurat_object = PEBBLES_soupx, group.by = "Samples", reduction = "pca")
ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/pca.pdf", 
       device = NULL,
       height = 8.5,
       width = 12)
ElbowPlot(PEBBLES_soupx)

#Filtered plots
Idents(PEBBLES_soupx) <- "Samples"
QC_Plot_UMIvsGene(seurat_object = PEBBLES_soupx, low_cutoff_gene = 600, high_cutoff_gene = 5000, low_cutoff_UMI = 500,
                  high_cutoff_UMI = 20000) +labs(x = "UMIs per nucleus", y = "Genes per nucleus" )
ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/postfilter_UMIvsGENE.pdf", 
       device = NULL,
       height = 8.5,
       width = 12)
QC_Plot_GenevsFeature(seurat_object = PEBBLES_soupx, feature1 = "percent.mito", low_cutoff_gene = 600,
                      high_cutoff_gene = 5000, high_cutoff_feature = 5)+labs(x = "Percent.mito per nucleus", y = "Genes per nucleus")
ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/postfilter_genevsmito.pdf", 
       device = NULL,
       height = 8.5,
       width = 12)

# Run KNN (K-nearest neighbor) followed by SNN (Shared-nearest neighbor)
PEBBLES_soupx <- FindNeighbors(PEBBLES_soupx, reduction="pca", dims = 1:20)
PEBBLES_soupx <- FindClusters(
  object = PEBBLES_soupx, 
  #resolution = seq(0.25,4,0.25), 
  resolution = 0.10,
  verbose = FALSE
)
sapply(grep("res",colnames(PEBBLES_soupx@meta.data),value = TRUE),
       function(x) length(unique(PEBBLES_soupx@meta.data[,x])))
Idents(PEBBLES_soupx) <- "RNA_snn_res.0.1"
table(PEBBLES_soupx$predicted.id,PEBBLES_soupx$RNA_snn_res.0.1)

#exclude the outlier sample
#PEBBLES_soupx <- test_object[, test_object@meta.data$Samples != "36_PCB_HET"]

#Cluster cells into UMAP
PEBBLES_soupx <- RunUMAP(
  object = PEBBLES_soupx,
  reduction.use = "pca",
  dims= 1:20,
  do.fast = TRUE,
  seed.use = 1234)
DimPlot_scCustom(seurat_object = PEBBLES_soupx, group.by = "Samples", reduction = "umap", pt.size = 0.5)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/UMAP_samples_regressed.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
DimPlot_scCustom(seurat_object = PEBBLES_soupx, group.by = "predicted.id", reduction = "umap", pt.size = 0.5)
Idents(PEBBLES_soupx) <- "predicted.id"
levels(PEBBLES_soupx) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

DimPlot_scCustom(seurat_object = PEBBLES_soupx, pt.size = 0.5, figure_plot = TRUE, label = FALSE)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/cell_type_UMAP.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
DimPlot_scCustom(seurat_object = PEBBLES_soupx, reduction = "umap", pt.size = 0.5, group.by = 'RNA_snn_res.0.25', figure_plot = TRUE)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/snn.res.0.25.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

Idents(PEBBLES_soupx) <- "RNA_snn_res.0.25"
PEBBLES_soupx <- BuildClusterTree(
  PEBBLES_soupx, dims = 1:20)
PlotClusterTree(PEBBLES_soupx)

#Find celltype markers for validation
all_markers <- FindAllMarkers(object = PEBBLES_soupx) %>%
  Add_Pct_Diff() %>%
  filter(pct_diff > 0.6)

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, named_vector = FALSE,
                                   make_unique = TRUE)

Clustered_DotPlot(seurat_object = PEBBLES_soupx, features = top_markers)
DotPlot_scCustom(seurat_object = PEBBLES_soupx, features = top_markers, flip_axes = T, remove_axis_titles = FALSE)
DoHeatmap(object = PEBBLES_soupx, features = top_markers, group.by = c("Treatment", "Samples"))
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/celltype_markers.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
Idents(PEBBLES_soupx) <- "Treatment"
DimPlot_scCustom(seurat_object = PEBBLES_soupx, pt.size = 0.5, figure_plot = FALSE, group.by = "Treatment")
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/Treatment.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
DimPlot_scCustom(seurat_object = PEBBLES_soupx, pt.size = 0.5, figure_plot = FALSE, group.by = "Genotype", )
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/Genotype.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
FeaturePlot_scCustom(PEBBLES_soupx, features = "Mecp2", pt.size = 0.05)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Rett_PCB_snRNAseq/Individual_plots/PCB_Mecp2.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
FeaturePlot_scCustom(PEBBLES_soupx, features = "Mecp2", pt.size = 0.05)

save(PEBBLES_soupx, file="/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_soupx_labeled.RData")
