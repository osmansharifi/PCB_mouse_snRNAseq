# load libraries
packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "dwtools", "devtools", "dplyr", "patchwork", "scCustomize")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

#Load data
setwd("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/")
load("../06_hdWGCNA/PEBBLES_clean.RData")
counts_table <- read.table("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/PEBBLES_Mecp2_alles.txt", sep="\t", header=FALSE)
#headers for the txt file
names(counts_table) <- c("Barcode", "UMI", "WT", "MUT", "Body", "Sample")
head(counts_table)

#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")


Mecp2_wt_mut_counts = counts_table
# Adding body counts into WT or MUT counts
new.counts.wt = sapply(1:nrow(Mecp2_wt_mut_counts), function(x) {
  ifelse(Mecp2_wt_mut_counts$Body[x]>0 & Mecp2_wt_mut_counts$WT[x]>0, 
         Mecp2_wt_mut_counts$Body[x]+Mecp2_wt_mut_counts$WT[x],
         Mecp2_wt_mut_counts$WT[x])
})

new.counts.mut = sapply(1:nrow(Mecp2_wt_mut_counts), function(x) {
  ifelse(Mecp2_wt_mut_counts$Body[x]>0 & Mecp2_wt_mut_counts$MUT[x]>0, 
         Mecp2_wt_mut_counts$Body[x]+Mecp2_wt_mut_counts$MUT[x],
         Mecp2_wt_mut_counts$MUT[x])
})

Mecp2_wt_mut_counts$WT = new.counts.wt
Mecp2_wt_mut_counts$MUT = new.counts.mut

# checking to make sure no duplicate barcodes - answer should be true
is.unique(Mecp2_wt_mut_counts$Barcode)

# making a data frame with just one vector representing all barcodes present in Seurat object by removing everything after the - symbol in the barcode names in the Seurat object
dataFrame = data.frame(Barcode = sub("-.*", "", names(Idents(PEBBLES_soupx))))

# adding in an id column to make sure we can sort the data back to the original order found in the Seurat Object
dataFrame$id = 1:nrow(dataFrame)

#checking to make sure all barcodes in Seurat Object are unique
is.unique(dataFrame$Barcode)

Mecp2_wt_mut_counts$Sample <- gsub("/$", "", Mecp2_wt_mut_counts$Sample)
Mecp2_wt_mut_counts$Barcode <- paste(Mecp2_wt_mut_counts$Sample, Mecp2_wt_mut_counts$Barcode, sep = "_")
is.unique(Mecp2_wt_mut_counts$Barcode)
# merging Mecp2 counts with Barcodes in the order that they appear in the Seurat Object
merged = merge(dataFrame, Mecp2_wt_mut_counts, by="Barcode", all=T)

# reordering to original Seurat Object order
merged = merged[order(merged$id),]

# converting NAs to zeros
merged[is.na(merged)] <- 0

# Number of barcodes that are not in Seurat object
number = nrow(merged) - PEBBLES_soupx@assays$RNA@counts@Dim[2] # x number of Barcodes in the Mecp2 WT vs. MUT counts file did not have corresponding Barcodes in the Seurat Object

# Removing rows with barcodes that are not in Seurat object
merged = merged[c(1:(nrow(merged)-number)),]

nrow(merged) # should equal length of all metadata 

### Adding counts to Seurat metadata 

PEBBLES_soupx$WT_Mecp2 = merged$WT

PEBBLES_soupx$MUT_Mecp2 = merged$MUT

### Adding counts to Seurat counts as two new genes

PEBBLES_soupx@assays$RNA@counts = rbind(PEBBLES_soupx@assays$RNA@counts, merged$WT, merged$MUT)
nrow(PEBBLES_soupx@assays$RNA@counts)
rownames(PEBBLES_soupx@assays$RNA@counts) = c(rownames(PEBBLES_soupx@assays$RNA@counts)[-c(18612:18613)], "Mecp2-WT", "Mecp2-MUT")

# checking to make sure they are added
PEBBLES_soupx@assays$RNA@counts[18612:18613, 1:5]

#E18 <- subset(x = PEBBLES_soupx, subset = orig.ident == c("MUT_F_E18_WB1", "MUT_F_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2"))
PEBBLES_soupx <- NormalizeData(PEBBLES_soupx)
PEBBLES_soupx <- ScaleData(PEBBLES_soupx)
PEBBLES_soupx <- FindVariableFeatures(PEBBLES_soupx, selection.method = "vst", nfeatures = 3000)
PEBBLES_soupx <- RunPCA(object = PEBBLES_soupx, verbose = FALSE)
PEBBLES_soupx <- RunUMAP(object = PEBBLES_soupx, dims = 1:20, verbose = FALSE)
PEBBLES_soupx <- FindNeighbors(object = PEBBLES_soupx, dims = 1:20, verbose = FALSE)
PEBBLES_soupx <- FindClusters(object = PEBBLES_soupx, verbose = FALSE)
DimPlot(object = PEBBLES_soupx, label = TRUE, group.by = 'celltype.call') + NoLegend() + ggtitle("sctransform")# saving new Seurat object
FeaturePlot_scCustom(seurat_object = PEBBLES_soupx, features = 'WT_Mecp2')
FeaturePlot_scCustom(seurat_object = PEBBLES_soupx, features = 'MUT_Mecp2')
FeaturePlot_scCustom(seurat_object = PEBBLES_soupx, features = 'Mecp2', split.by = "Sex")

#Function counts the percent of total cells that express specific genes
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

PrctCellExpringGene(PEBBLES_soupx, genes =c("WT-Mecp2","MUT-Mecp2"), group.by = "all")

GOI1 <- 'WT_Mecp2' #you will have to name your first gene here, im choosing PDX1 as an example
GOI2 <- 'MUT_Mecp2' #you will have to name your second gene here, im choosing INS as an example
GOI1.cutoff <- 1 #Assumption: gene count cutoff is 1, assuming atleast 1 count is REAL
GOI2.cutoff <- 1 #Assumption: gene count cutoff is 1, assuming atleast 1 count is REAL

# Time to party
GOI1.cells <- length(which(FetchData(E18, vars = GOI1) > GOI1.cutoff))
GOI2.cells <- length(which(FetchData(E18, vars = GOI2) > GOI2.cutoff))
GOI1_GOI2.cells <- length(which(FetchData(E18, vars = GOI2) > GOI2.cutoff & FetchData(E18, vars = GOI1) > GOI1.cutoff))
all.cells.incluster <- table(E18$orig.ident)
GOI1.cells/all.cells.incluster*100 # Percentage of cells in Beta that express GOI1
GOI2.cells/all.cells.incluster*100 #Percentage of cells in Beta that express GOI2
GOI1_GOI2.cells/all.cells.incluster*100 #Percentage of cells in Beta that co-express GOI1 + GOI2

#Marker Genes
cell_markers_manual <- c("Plch2","Sst","Vip", "Pvalb", "Slc17a8", "Macc1", "Rorb", "Fezf2", "Rprm", "Aqp4", "Rassf10", "Kcnj8", "Slc17a7", "Gad2", "Aspa")

FeaturePlot_scCustom(seurat_object = p120, features = cell_markers_manual)
DotPlot_scCustom(seurat_object = p120, features = cell_markers_manual, group.by = "celltype.call", x_lab_rotate = TRUE, flip_axes = TRUE)
ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/celltype_markers.pdf")
Clustered_DotPlot(seurat_object = p120, features = cell_markers_manual, colors_use_idents = TRUE, row_label_size = 14)
ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/celltype_markers_clustered.pdf")


markers_df <- FindAllMarkers(
  object = p120, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  thresh.use = 0.25
)
markers_all_single <- markers_df[markers_df$gene %in% names(table(markers_df$gene))[table(markers_df$gene) == 1],]
top5 <- markers_all_single %>% group_by(cluster) %>% top_n(5, avg_log2FC)
dim(top5)
DoHeatmap(
  object = p120, 
  features = top5$gene,
  labels = FALSE) 

save(p120, file="/Users/osman/Desktop/LaSalle_lab/Seurat_objects/p120.RData")
save(p120, file=glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Parsing_Mecp2_trnx_expression/{p120.name}/{p120.name}_with_Mecp2_WT_MUT.RData"))

