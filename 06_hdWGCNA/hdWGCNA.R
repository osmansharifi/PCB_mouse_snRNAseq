# This program will perform hdWGCNA analysis on single cell data
# single-cell analysis package
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(dplyr)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(1234)

# load the snRNA-seq dataset
load("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_soupx_labeled.RData")
setwd("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/06_hdWGCNA")
PEBBLES_soupx$cell_type = PEBBLES_soupx$predicted.id

# Select columns that do not start with 'predicted' or start with 'predicted.id'
# Find columns starting with 'prediction' or 'RNA'
cols_to_null <- grep("^prediction|^RNA", names(PEBBLES_soupx@meta.data), value = TRUE)

# Set those columns to NULL
PEBBLES_soupx@meta.data[, cols_to_null] <- NULL

# Set up multithreading
enableWGCNAThreads(nThreads = 16)
allowWGCNAThreads(nThreads = 16)

# Prepare Seurat Object for WGCNA
metadata <- PEBBLES_soupx@meta.data
timepoint <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[3])
})
PEBBLES_soupx@meta.data$Time_Point <- unlist(timepoint)

genotype <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[1])
})
PEBBLES_soupx@meta.data$genotype <- unlist(genotype)


# Preprocess
PEBBLES_soupx <- NormalizeData(PEBBLES_soupx, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(PEBBLES_soupx)
PEBBLES_soupx <- ScaleData(PEBBLES_soupx)#, features = all.genes)
PEBBLES_soupx <- FindVariableFeatures(PEBBLES_soupx, selection.method = "vst", nfeatures = 2000)
PEBBLES_soupx <- RunPCA(PEBBLES_soupx, features = VariableFeatures(object = PEBBLES_soupx))
PEBBLES_soupx <- RunUMAP(PEBBLES_soupx, dims = 1:20)
DimPlot(PEBBLES_soupx, group.by='celltype.call', label=TRUE) +
  umap_theme() 
DefaultAssay(PEBBLES_soupx) <- 'RNA'
Idents(PEBBLES_soupx) <- 'celltype.call'
# Set up Seurat object for WGCNA
PEBBLES_soupx <- SetupForWGCNA(
  PEBBLES_soupx,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.1, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "postnatal_mouse_cortex" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
PEBBLES_soupx <- MetacellsByGroups(
  PEBBLES_soupx,
  group.by = c("cell_type", "Treatment"), # specify the columns in PEBBLES_soupx@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
PEBBLES_soupx <- NormalizeMetacells(seurat_obj = PEBBLES_soupx, wgcna_name = "postnatal_mouse_cortex",)
Idents(PEBBLES_soupx) <- "cell_type"
# Set up the expression matrix
PEBBLES_soupx <- SetDatExpr(
  PEBBLES_soupx,
  group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Non-neuronal", "Astro", "Oligo", "Lamp5"), # the name of the group of interest in the group.by column
  group.by="cell_type", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
PEBBLES_soupx <- TestSoftPowers(
  PEBBLES_soupx,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(PEBBLES_soupx)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
ggplot2::ggsave("Softpowerthreshold.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# construct co-expression network:
PEBBLES_soupx <- ConstructNetwork(
  PEBBLES_soupx, soft_power=8,
  setDatExpr=FALSE,
  overwrite_tom = TRUE# name of the topoligical overlap matrix written to disk
)

PlotDendrogram(PEBBLES_soupx, main='hdWGCNA Dendrogram')
ggplot2::ggsave("WGCNA_Dendrogram.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# need to run ScaleData first or else harmony throws an error:
PEBBLES_soupx <- ScaleData(PEBBLES_soupx, features=VariableFeatures(PEBBLES_soupx))

# compute all MEs in the full single-cell dataset
PEBBLES_soupx <- ModuleEigengenes(
  PEBBLES_soupx,
  group.by.vars="orig.ident"
)

# harmonized module eigengenes:
hMEs <- GetMEs(PEBBLES_soupx)

# module eigengenes:
MEs <- GetMEs(PEBBLES_soupx, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
PEBBLES_soupx <- ModuleConnectivity(
  PEBBLES_soupx,
  group.by = 'celltype.call', group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Sncg", "Non-neuronal", "Astro", "Oligo", "Lamp5")
)

# plot genes ranked by kME for each module
p <- PlotKMEs(PEBBLES_soupx, ncol=5)

p
ggplot2::ggsave("kME.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
PEBBLES_soupx <- ModuleExprScore(
  PEBBLES_soupx,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
PEBBLES_soupx <- ModuleExprScore(
  PEBBLES_soupx,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  PEBBLES_soupx,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)
ggplot2::ggsave("module_UMAPs.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

levels(PEBBLES_soupx) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
DimPlot_scCustom(seurat_object = PEBBLES_soupx, label = FALSE, pt.size = 0.5)
ggplot2::ggsave("celltype_UMAPs.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  PEBBLES_soupx,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE) # depending on Seurat vs UCell for gene scoring


# stitch together with patchwork
wrap_plots(plot_list, ncol=3)
ggplot2::ggsave("hubgene_scores_UMAPs.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# plot module correlagram
ModuleCorrelogram(PEBBLES_soupx)
ggplot2::ggsave("module_to_module_cor.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

# get hMEs from seurat object
MEs <- GetMEs(PEBBLES_soupx, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
PEBBLES_soupx@meta.data <- cbind(PEBBLES_soupx@meta.data, MEs)
# plot with Seurat's DotPlot function
DotPlot_scCustom(seurat_object = PEBBLES_soupx, features=mods, flip_axes = TRUE, x_lab_rotate = TRUE, remove_axis_titles = FALSE) + xlab("Modules") + ylab("Cell_Type")
ggplot2::ggsave("Average_expression_hubgenes.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

## Compute Correlations
# convert genotype to factor
PEBBLES_soupx$Genotype <- as.factor(PEBBLES_soupx$genotype)
# convert time point to factor
PEBBLES_soupx$Time_Point <- as.factor(PEBBLES_soupx$Time_Point)
# convert celltype to factor
PEBBLES_soupx$celltype.call <- as.factor(PEBBLES_soupx$celltype.call)
# convert sample name to factor
PEBBLES_soupx$orig.ident <- as.factor(PEBBLES_soupx$orig.ident)
# convert sex to factor
PEBBLES_soupx$Sex <- as.factor(PEBBLES_soupx$Sex)
# convert sex to factor
PEBBLES_soupx$Disease_score <- as.numeric(PEBBLES_soupx$Disease_score)
# convert bodyweight to factor
PEBBLES_soupx$Body_weight <- as.numeric(PEBBLES_soupx$Body_weight)

# list of traits to correlate
cur_traits <- c('Sex','Time_Point','Genotype', 'Disease_score', 'Body_weight')

PEBBLES_soupx <- ModuleTraitCorrelation(
  PEBBLES_soupx,
  traits = cur_traits,
  group.by='celltype.call'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(PEBBLES_soupx)

names(mt_cor$cor)
PlotModuleTraitCorrelation(
  PEBBLES_soupx,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 4,
  text_digits = 4,
  text_color = 'black',
  high_color = '#B2182B',
  mid_color = '#EEEEEE',
  low_color = '#2166AC',
  plot_max = 0.4,
  combine=TRUE
)

#Warning messages:
#1: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits,  :Trait Sex is a factor with levels Female, Male. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#2: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits,  :Trait Time_Point is a factor with levels P30, P60, P120, P150. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#3: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits,  :Trait Genotype is a factor with levels MUT, WT. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#4: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits,  :Trait disease_score is a factor with levels 0, 0.5, 1, 2, 3.5, 5. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#5: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits,  :Trait body_weight is a factor with levels 22, 24, 25, 29, 30, 35, 38, 45, 47, 50. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?

mod_trait_cor <- mt_cor$cor
                                                                                           
# get modules
modules <- GetModules(PEBBLES_soupx)
head(modules)
write.csv(modules, "modules.csv", row.names = FALSE)
# get hub genes
hub_genesdf <- GetHubGenes(PEBBLES_soupx, n_hubs = 10)
head(hub_genesdf)
write.csv(hub_genesdf, "top10_hub_genes.csv", row.names = FALSE)

save(list = ls(), file = "mouse_cortex_WGCNA.RData")
