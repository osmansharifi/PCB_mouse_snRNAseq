# This program will perform hdWGCNA analysis on single cell data
# single-cell analysis package
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(UCell)
library(magrittr)
library(igraph)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(1234)

# load the snRNA-seq dataset
setwd("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/06_hdWGCNA")
load("PEBBLES_clean.RData")

# Set up multithreading
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
DefaultAssay(PEBBLES_soupx) <- 'RNA'
Idents(PEBBLES_soupx) <- 'cell_type'
# Set up Seurat object for WGCNA
PEBBLES_soupx <- SetupForWGCNA(
  PEBBLES_soupx,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.15, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "pebbles_cortex_hdwgcna" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
PEBBLES_soupx <- MetacellsByGroups(
  PEBBLES_soupx,
  group.by = c("cell_type"), # specify the columns in PEBBLES_soupx@meta.data to group by
  wgcna_name = 'pebbles_cortex_hdwgcna', 
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
PEBBLES_soupx <- NormalizeMetacells(seurat_obj = PEBBLES_soupx, wgcna_name = "pebbles_cortex_hdwgcna",)
Idents(PEBBLES_soupx) <- "cell_type"
# Set up the expression matrix
PEBBLES_soupx <- SetDatExpr(
  PEBBLES_soupx,
  group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Non-neuronal", "Astro", "Oligo", "Lamp5", "Sncg"), # the name of the group of interest in the group.by column
  group.by="cell_type", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
PEBBLES_soupx <- TestSoftPowers(
  PEBBLES_soupx,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
) # this errors out if there are columns that are constant and dont change 

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
  PEBBLES_soupx, soft_power=5,
  setDatExpr=FALSE,
  overwrite_tom = TRUE# name of the topoligical overlap matrix written to disk
)

PlotDendrogram(PEBBLES_soupx, main='hdWGCNA PEBBLES Dendrogram')
ggplot2::ggsave("WGCNA_Dendrogram.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

#PEBBLES_soupx@misc$pebbles_cortex_hdwgcna$wgcna_modules$module <- paste0(PEBBLES_soupx@misc$pebbles_cortex_hdwgcna$wgcna_modules$module, "_")
# need to run ScaleData first or else harmony throws an error:
PEBBLES_soupx <- ScaleData(PEBBLES_soupx, features=VariableFeatures(PEBBLES_soupx))

# compute all MEs in the full single-cell dataset
PEBBLES_soupx <- ModuleEigengenes(
  PEBBLES_soupx,
  group.by.vars="Group",
  wgcna_name = "pebbles_cortex_hdwgcna",
  verbose = TRUE,
  pc_dim = c(1:20)
)
PEBBLES_soupx <- ModuleEigengenes(
  PEBBLES_soupx)
# harmonized module eigengenes:
hMEs <- GetMEs(PEBBLES_soupx)

# module eigengenes:
MEs <- GetMEs(PEBBLES_soupx, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
PEBBLES_soupx <- ModuleConnectivity(
  PEBBLES_soupx,
  group.by = 'cell_type', group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Sncg", "Non-neuronal", "Astro", "Oligo", "Lamp5")
)

# plot genes ranked by kME for each module
PlotKMEs(PEBBLES_soupx, ncol=3)
ggplot2::ggsave("PEBBLES_kME.pdf",
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
ggplot2::ggsave("PEBBLES_module_UMAPs.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

levels(PEBBLES_soupx) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
DimPlot_scCustom(seurat_object = PEBBLES_soupx, label = FALSE, pt.size = 0.5, figure_plot = TRUE)
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
ggplot2::ggsave("PEBBLES_hubgene_scores_UMAPs.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# plot module correlagram
ModuleCorrelogram(PEBBLES_soupx)
ggplot2::ggsave("PEBBLES_module_to_module_cor.pdf",
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
ggplot2::ggsave("PEBBLES_Average_expression_hubgenes.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

################################
## Add Traits to the metadata ##
################################
traits_sheet <- read.table("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/06_hdWGCNA/PEBBLES_traits.csv", sep = ",", header = TRUE)
traits_sheet$X.1 <- NULL
traits_sheet$X <- NULL
traits_sheet$X.2 <- NULL
traits_sheet$X.3 <- NULL
traits_sheet$X.4 <- NULL
traits_sheet$X.5 <- NULL

# Get the sample names from the Seurat object
sample_names <- PEBBLES_soupx@meta.data$Samples

# Add the exposure, weight and pregnant columns from traits_sheet
PEBBLES_soupx@meta.data$Exposure_duration <- NA
PEBBLES_soupx@meta.data$Weight <- NA
PEBBLES_soupx@meta.data$Pregnant <- NA
# Match samples and fill in the data
for (i in 1:length(sample_names)) {
  sample_name <- sample_names[i]
  row_index <- traits_sheet$Samples == sample_name
  if (sum(row_index) != 0) {
    PEBBLES_soupx@meta.data$Exposure_duration[i] <- traits_sheet$Exposure_duration[row_index]
    PEBBLES_soupx@meta.data$Weight[i] <- traits_sheet$Weight[row_index]
    PEBBLES_soupx@meta.data$Pregnant[i] <- traits_sheet$Pregnant[row_index]
  }
}

##########################
## Compute Correlations ##
##########################
# set as factor or numeric
PEBBLES_soupx$Genotype <- as.factor(PEBBLES_soupx$Genotype)
PEBBLES_soupx$Treatment <- as.factor(PEBBLES_soupx$Treatment)
PEBBLES_soupx$Exposure_duration <- as.numeric(PEBBLES_soupx$Exposure_duration)
PEBBLES_soupx$Weight <- as.numeric(PEBBLES_soupx$Weight)
PEBBLES_soupx$Pregnant <- as.factor(PEBBLES_soupx$Pregnant)
PEBBLES_soupx$Samples <- as.factor(PEBBLES_soupx$Samples)
PEBBLES_soupx$Group <- paste(PEBBLES_soupx$Genotype, PEBBLES_soupx$Treatment, sep = "-")
PEBBLES_soupx$Group <- as.factor(PEBBLES_soupx$Group)
# list of traits to correlate
cur_traits <- c('Genotype', 'Treatment', 'Exposure_duration', 'Weight', 'Pregnant', 'Group')

PEBBLES_soupx <- ModuleTraitCorrelation(
  PEBBLES_soupx,
  traits = cur_traits,
  group.by='cell_type'
)

#Warning messages:
#  1: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits, group.by = "cell_type") :
#  Trait Samples is a factor with levels 24_PCB_WT, 25_VEHICLE_WT, 27_PCB_HET, 27_PCB_HET_2, 28_VEHICLE_HET, 29_VEHICLE_WT, 30_VEHICLE_WT, #30_VEHICLE_WT_2, 31_PCB_WT, 37_PCB_WT, 37_PCB_WT_2, 38_VEHICLE_HET, 39_PCB_HET, 40_VEHICLE_HET, 40_VEHICLE_HET_2. Levels will be converted to #numeric IN THIS ORDER for the correlation, is this the expected order?
#  2: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits, group.by = "cell_type") :
#  Trait Genotype is a factor with levels HET, WT. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected #order?
#  3: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits, group.by = "cell_type") :
#  Trait Treatment is a factor with levels PCB, VEHICLE. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the #expected order?
#  4: In ModuleTraitCorrelation(PEBBLES_soupx, traits = cur_traits, group.by = "cell_type") :
#  Trait Pregnant is a factor with levels No, Yes. Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected #order?
# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(PEBBLES_soupx)

names(mt_cor$cor)
PlotModuleTraitCorrelation(
  PEBBLES_soupx,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 3,
  text_digits = 4,
  text_color = 'black',
  high_color = '#B2182B',
  mid_color = '#EEEEEE',
  low_color = '#2166AC',
  plot_max = 0.8,
  combine=TRUE
)

mod_trait_cor <- mt_cor$cor
                                                                                           
# get modules
modules <- GetModules(PEBBLES_soupx)
head(modules)
write.csv(modules, "modules.csv", row.names = FALSE)
# get hub genes
hub_genesdf <- GetHubGenes(PEBBLES_soupx, n_hubs = 10)
head(hub_genesdf)
write.csv(hub_genesdf, "top10_hub_genes.csv", row.names = FALSE)

save(list = ls(), file = "PEBBLES_cortex_WGCNA.RData")

ModuleNetworkPlot(
  PEBBLES_soupx,
  outdir = 'ModuleNetworks'
)

# hubgene network
HubGeneNetworkPlot(
  PEBBLES_soupx,
  n_hubs = 1, n_other=149,
  edge_prop = 1,
  mods = 'all',
  edge.alpha = 0.5,
  vertex.label.cex = 1,
  hub.vertex.size = 6
)

g <- HubGeneNetworkPlot(PEBBLES_soupx,  return_graph=TRUE)

seurat_obj <- RunModuleUMAP(
  PEBBLES_soupx,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(PEBBLES_soupx)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  PEBBLES_soupx,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)

# Add a column for module names
test <- cbind(mod_trait_cor$all_cells, mod_trait_cor$Astro)

dbs <- "KEGG_2019_Mouse"
# perform enrichment tests
PEBBLES_soupx <- RunEnrichr(
  PEBBLES_soupx,
  dbs="KEGG_2019_Mouse", # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(PEBBLES_soupx)

# make GO term plots:
EnrichrBarPlot(
  PEBBLES_soupx,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
EnrichrDotPlot(
  PEBBLES_soupx,
  mods = "all", # use all modules (this is the default behavior)
  database = dbs, # this has to be one of the lists we used above!!!
  n_terms=5 # number of terms for each module
)
ggplot2::ggsave("PEBBLES_top5_KEGG.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

# compute cell-type marker genes with Seurat:
Idents(PEBBLES_soupx) <- PEBBLES_soupx$cell_type
markers <- Seurat::FindAllMarkers(
  PEBBLES_soupx,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  PEBBLES_soupx,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)
# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
wrap_plots(plot_list, ncol=3)
ggplot2::ggsave("PEBBLES_module_odds.ratio.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# plot odds ratio of the overlap as a dot plot
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
ggplot2::ggsave("PEBBLES_cellmarker_module_overlap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

# Extract the relevant information from meta.data
library(ggpubr)
plot_data <- PEBBLES_soupx@meta.data[, c("Group", "brown", "Exposure_duration", "Weight", "Pregnant")]
plot_data$Group <- factor(plot_data$Group, levels = c("WT-VEHICLE", "WT-PCB", "HET-VEHICLE", "HET-PCB"))
p <- ggplot(plot_data, aes(x = Group, y = brown, fill = Pregnant)) +
  geom_violin() +
  labs(
    title = "Violin Plot of Group vs brownE",
    x = "Group",
    y = "Module Eigengene"
  ) +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'Module Trait correlation') +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 90, size = 12, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 8, face = 'plain', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),
    # removes the border
    legend.key.size = unit(1, "cm"),
    # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"),
    # Text size
    title = element_text(size = 18, face = "bold")
  )
p + stat_compare_means(
  comparisons = list(c("HET-PCB", "WT-VEHICLE"), c("HET-VEHICLE", "WT-VEHICLE"), c("WT-PCB", "WT-VEHICLE")),
  method = "t.test",
  label = "p.format"
)
ggplot2::ggsave("PEBBLES_module_pregnant_corr.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

p <- ggplot(plot_data, aes(x = Group, y = brown, fill = Exposure_duration)) +
  geom_violin() +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'brown Module Trait correlation', y = "Module Eigengene") +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 90, size = 12, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 8, face = 'plain', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),
    # removes the border
    legend.key.size = unit(1, "cm"),
    # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"),
    # Text size
    title = element_text(size = 18, face = "bold")
  )
p + stat_compare_means(
  comparisons = list(c("HET-PCB", "WT-VEHICLE"), c("HET-VEHICLE", "WT-VEHICLE"), c("WT-PCB", "WT-VEHICLE")),
  method = "t.test",
  label = "p.format"
)
ggplot2::ggsave("PEBBLES_module_ExposureDuration_corr.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

p <- ggplot(plot_data, aes(x = Group, y = brown, fill = Pregnant)) +
  geom_violin() +
  labs(
    title = "Violin Plot of Group vs brownE",
    x = "Group",
    y = "Module Eigengene"
  ) +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'brown Module Trait correlation', y = "Module Eigengene") +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 90, size = 12, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 8, face = 'plain', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),
    # removes the border
    legend.key.size = unit(1, "cm"),
    # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"),
    # Text size
    title = element_text(size = 18, face = "bold")
  )
p + stat_compare_means(
  comparisons = list(c("HET-PCB", "WT-VEHICLE"), c("HET-VEHICLE", "WT-VEHICLE"), c("WT-PCB", "WT-VEHICLE")),
  method = "t.test",
  label = "p.format"
)
ggplot2::ggsave("PEBBLES_module_pregnant_corr.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

brown_cor <- data.frame(mt_cor$cor$Sst[,5:5])
brown_fdr <- data.frame(mt_cor$fdr$Sst[,5:5])
brown <- merge(brown_cor, brown_fdr, by = "row.names", all = TRUE)


