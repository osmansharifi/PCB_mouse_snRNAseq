# This program will perform hdWGCNA analysis on single cell data
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(1234)

# load the snRNA-seq dataset
load("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_soupx_labeled.RData")
setwd("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis")
Idents(all.cortex.combined) <- "celltype.call"
# Set up multithreading
enableWGCNAThreads(nThreads = 16)

adult_postnatal <- subset(all.cortex.combined, subset = Age != 'E18')
# Prepare Seurat Object for WGCNA
metadata <- adult_postnatal@meta.data
timepoint <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[3])
})
adult_postnatal@meta.data$Time_Point <- unlist(timepoint)

genotype <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[1])
})
adult_postnatal@meta.data$genotype <- unlist(genotype)


# Preprocess
adult_postnatal <- NormalizeData(adult_postnatal, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(adult_postnatal)
adult_postnatal <- ScaleData(adult_postnatal)#, features = all.genes)
adult_postnatal <- FindVariableFeatures(adult_postnatal, selection.method = "vst", nfeatures = 2000)
adult_postnatal <- RunPCA(adult_postnatal, features = VariableFeatures(object = adult_postnatal))
adult_postnatal <- RunUMAP(adult_postnatal, dims = 1:20)
DimPlot(adult_postnatal, group.by='celltype.call', label=TRUE) +
  umap_theme() 
DefaultAssay(adult_postnatal) <- 'RNA'
Idents(adult_postnatal) <- 'celltype.call'
# Set up Seurat object for WGCNA
adult_postnatal <- SetupForWGCNA(
  adult_postnatal,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "postnatal_mouse_cortex" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
adult_postnatal <- MetacellsByGroups(
  adult_postnatal,
  group.by = c("celltype.call", "Time_Point", "Sex", "genotype"), # specify the columns in adult_postnatal@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype.call' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
adult_postnatal <- NormalizeMetacells(adult_postnatal)
Idents(adult_postnatal) <- "celltype.call"
# Set up the expression matrix
adult_postnatal <- SetDatExpr(
  adult_postnatal,
  group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Sncg", "Non-neuronal", "Astro", "Oligo", "Lamp5"), # the name of the group of interest in the group.by column
  group.by="celltype.call", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
adult_postnatal <- TestSoftPowers(
  adult_postnatal,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(adult_postnatal)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
ggplot2::ggsave("Softpowerthreshold.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# construct co-expression network:
adult_postnatal <- ConstructNetwork(
  adult_postnatal, soft_power=8,
  setDatExpr=FALSE,
  overwrite_tom = TRUE# name of the topoligical overlap matrix written to disk
)

PlotDendrogram(adult_postnatal, main='hdWGCNA Dendrogram')
ggplot2::ggsave("WGCNA_Dendrogram.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# need to run ScaleData first or else harmony throws an error:
adult_postnatal <- ScaleData(adult_postnatal, features=VariableFeatures(adult_postnatal))

# compute all MEs in the full single-cell dataset
adult_postnatal <- ModuleEigengenes(
  adult_postnatal,
  group.by.vars="orig.ident"
)

# harmonized module eigengenes:
hMEs <- GetMEs(adult_postnatal)

# module eigengenes:
MEs <- GetMEs(adult_postnatal, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
adult_postnatal <- ModuleConnectivity(
  adult_postnatal,
  group.by = 'celltype.call', group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Sncg", "Non-neuronal", "Astro", "Oligo", "Lamp5")
)

# plot genes ranked by kME for each module
p <- PlotKMEs(adult_postnatal, ncol=5)

p
ggplot2::ggsave("kME.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
adult_postnatal <- ModuleExprScore(
  adult_postnatal,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
adult_postnatal <- ModuleExprScore(
  adult_postnatal,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  adult_postnatal,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)
ggplot2::ggsave("module_UMAPs.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

levels(adult_postnatal) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
DimPlot_scCustom(seurat_object = adult_postnatal, label = FALSE, pt.size = 0.5)
ggplot2::ggsave("celltype_UMAPs.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  adult_postnatal,
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
ModuleCorrelogram(adult_postnatal)
ggplot2::ggsave("module_to_module_cor.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

# get hMEs from seurat object
MEs <- GetMEs(adult_postnatal, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
adult_postnatal@meta.data <- cbind(adult_postnatal@meta.data, MEs)
# plot with Seurat's DotPlot function
DotPlot_scCustom(seurat_object = adult_postnatal, features=mods, flip_axes = TRUE, x_lab_rotate = TRUE, remove_axis_titles = FALSE) + xlab("Modules") + ylab("Cell_Type")
ggplot2::ggsave("Average_expression_hubgenes.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

## Compute Correlations
# convert genotype to factor
adult_postnatal$Genotype <- as.factor(adult_postnatal$genotype)
# convert time point to factor
adult_postnatal$Time_Point <- as.factor(adult_postnatal$Time_Point)
# convert celltype to factor
adult_postnatal$celltype.call <- as.factor(adult_postnatal$celltype.call)
# convert sample name to factor
adult_postnatal$orig.ident <- as.factor(adult_postnatal$orig.ident)
# convert sex to factor
adult_postnatal$Sex <- as.factor(adult_postnatal$Sex)
# convert sex to factor
adult_postnatal$Disease_score <- as.numeric(adult_postnatal$Disease_score)
# convert bodyweight to factor
adult_postnatal$Body_weight <- as.numeric(adult_postnatal$Body_weight)

# list of traits to correlate
cur_traits <- c('Sex','Time_Point','Genotype', 'Disease_score', 'Body_weight')

adult_postnatal <- ModuleTraitCorrelation(
  adult_postnatal,
  traits = cur_traits,
  group.by='celltype.call'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(adult_postnatal)

names(mt_cor$cor)
PlotModuleTraitCorrelation(
  adult_postnatal,
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
#1: In ModuleTraitCorrelation(adult_postnatal, traits = cur_traits,  :Trait Sex is a factor with levels Female, Male. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#2: In ModuleTraitCorrelation(adult_postnatal, traits = cur_traits,  :Trait Time_Point is a factor with levels P30, P60, P120, P150. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#3: In ModuleTraitCorrelation(adult_postnatal, traits = cur_traits,  :Trait Genotype is a factor with levels MUT, WT. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#4: In ModuleTraitCorrelation(adult_postnatal, traits = cur_traits,  :Trait disease_score is a factor with levels 0, 0.5, 1, 2, 3.5, 5. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?
#5: In ModuleTraitCorrelation(adult_postnatal, traits = cur_traits,  :Trait body_weight is a factor with levels 22, 24, 25, 29, 30, 35, 38, 45, 47, 50. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?

mod_trait_cor <- mt_cor$cor
                                                                                           
# get modules
modules <- GetModules(adult_postnatal)
head(modules)
write.csv(modules, "modules.csv", row.names = FALSE)
# get hub genes
hub_genesdf <- GetHubGenes(adult_postnatal, n_hubs = 10)
head(hub_genesdf)
write.csv(hub_genesdf, "top10_hub_genes.csv", row.names = FALSE)

save(list = ls(), file = "mouse_cortex_WGCNA.RData")
