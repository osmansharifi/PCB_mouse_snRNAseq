# Tutorial: Deconvolution of mouse PatchSeq

This tutorial provides a guided deconvolution of the gene expression from mouse neurons profiled with PatchSeq produced by the Allen Brain Institute. Deconvolution requires a reference atlas to define the gene expression patterns of individual cell types. In this tutorial we will be using a matched mouse scRNAseq dataset available through the  [Allen Brain Atlas: Cell Types data portal](https://celltypes.brain-map.org/).

If you need to setup a python virtual enviroment refer to: [Python virtual env. tutorial](''). This tutorial must be run on `pggb` or another GPU node to utilize Tensorflow.


## Deconvolution goals
The following is a walkthrough of a standard deconvolution problem with `<NAME>` and has been designed to provide an overview of data preprocessing, deconvolution and sanity checking. Here, our primary goals include:

1. Setting up the two datasets and mapping marker genes from cluster names.
2. Train `<NAME>` on a large  reference atlas and performing deconvolution of mouse patchseq.
3. Perform standard sanity checks on estimated proportions and purified expression profiles.

## Data setup
The preprocessed patchseq and scRNAseq data for this tutorial can be obtained by running the R script: `/share/quonlab/wkdir/njjohans/public/deconv_tutorial/mouse_patchseq/mouse_scRNA_patchseq.R `.

Once you've generated the initial preprocessed data we can load it in to R:
```R
library(reticulate) ## Let's us call Python from R!
library(Seurat)
library(ComplexHeatmap)
library(circlize)
# library(patchSeqQC)
deconv = import("deconv") ## <NAME>!
options(stringsAsFactors=F)

## Load in helper functions/scripts
source('/share/quonlab/wkdir/njjohans/public/deconv_tutorial/scripts/nnDeconvClass.R')
source('/share/quonlab/wkdir/njjohans/public/deconv_tutorial/scripts/proportionHeatmap.R')
source('/share/quonlab/wkdir/njjohans/public/deconv_tutorial/scripts/proportionSanityCheck.R')

## User paths
working.dir = "/share/quonlab/wkdir/njjohans/public/deconv_tutorial/mouse_patchseq/" #where the data files are located

## Set the working directory
setwd(working.dir)

## Load in mouse (atlas) and mouse.patchseq Seurat objects.
load(paste0(working.dir, "mouse.rda"))

## Check namespace
ls()
```

### Cell annotation setup
Extract cell type annotations and group some cell types into higher level annotations.
```R
## Take a look at one of the cell annotations (there are a few levels of cell annotations in mouse@meta.data), we will focus on the higher level annotations for now.
print(table(mouse@meta.data$subclass_label))

## Ensure that cell type names follow standard python conventions, i.e. no "/" and " ".
celltype.deconv = gsub("/", "_", gsub(" ", "_", mouse@meta.data$subclass_label))

## Merge cell types to form high level annotations
celltype.deconv[which(celltype.deconv %in% c("CR", "Macrophage", "Meis2", "NP", "Serpinf1", "SMC", "VLMC"))] = "Non-neuronal"
celltype.deconv[grep("L5", celltype.deconv)] = "L5"
celltype.deconv[grep("L6", celltype.deconv)] = "L6"

## Store our updated cell type annotations in the Seurat object
mouse@meta.data$celltype.deconv = celltype.deconv
```
## scRNAseq atlas sub-sampling
Due to the size of the mouse scRNAseq atlas we need to sub-sample the data so that `<NAME>` will fit into GPU VRAM. We will subsample the cell types to contain at most 1000 cells: *(Note: This is not a standard step and will be removed after `<NAME>` has been updated to handle high memory deconvolution tasks.)*

```R
# Function to subsample cells
subsampleCells <- function(clustersF, subSamp=100, seed=5){
  # Returns a vector of TRUE false for choosing a maximum of subsamp cells in each cluster
  # clustersF = vector of cluster labels in factor format
  kpSamp = rep(FALSE,length(clustersF))
  for (cli in unique(as.character(clustersF))){
    set.seed(seed)
    seed   = seed+1
    kp     = which(clustersF==cli)
    kpSamp[kp[sample(1:length(kp),min(length(kp),subSamp))]] = TRUE
  }
  return(kpSamp)
}

## Keep 1000 random samples from each of our pre-defined celltype annotations
cells.keep = subsampleCells(mouse@meta.data$celltype.deconv, subSamp=1000)

## Subset the scRNAseq Seurat object
mouse = subset(mouse, cells=colnames(mouse)[cells.keep])
```

## Marker gene extraction
Next we determine a set of marker genes which maximally distinguish the cell types of interest by using the cluster annotations provided by the Allen Brain Institute.
```R
## dplyr makes extracting data from data.frames much easier
library(dplyr)

## Identify the genes in common between datasets
common.gene.set = intersect(rownames(mouse), rownames(mouse.patchseq))

## Celltype marker annotations
cluster.meta = data.frame(celltype=mouse@meta.data$celltype.deconv, marker=gsub(" ", "", gsub(" ", "_", mouse@meta.data$cluster_label)))

## Merge high level clustering using dplyr. *You may want to separate this and run line by line.*
cluster.markers = cluster.meta %>%
                    group_by(celltype) %>%
                    dplyr::summarise(markers=paste(unique(unlist(strsplit(marker, "_"))), collapse=" ")) %>%
                    mutate(marker.genes=strsplit(markers," "))

## Extract marker genes embedded in cell type annotation `cluster_label`
cluster.marker.list = lapply(as.list(cluster.markers$markers), function(x) unlist(strsplit(x, " ")))
cluster.marker.list = lapply(cluster.marker.list, function(x) x[which(x %in% common.gene.set)])
names(cluster.marker.list) = cluster.markers$celltype

## Finally we flatten the list into a data.frame
marker.anno.df = data.frame(celltype=character(), gene=character())
for(celltype in names(cluster.marker.list)){
  for(gene in cluster.marker.list[[celltype]]){
      marker.anno.df[nrow(marker.anno.df)+1,] = c(celltype, gene)
  }
}

## Get the unique markers, some clusters can share marker genes!
marker.genes = unique(marker.anno.df$gene)
```

## Atlas exploration
Let's take a quick look at our atlas now that we have marker genes and cell type calls!
```R
## Run dimensionality reduction on just the marker genes
mouse <- RunPCA(mouse, npcs = 30, features=marker.genes, verbose = FALSE)
mouse <- RunUMAP(mouse, reduction = "pca", dims = 1:30)

## Plot
umap.plot <- DimPlot(mouse, reduction = "umap", group.by = "celltype.deconv", label=TRUE, repel=TRUE)

## Save our plot
png("./mouse_scRNA_atlas.png", width=12, height=12, res=150,  units="in")
plot(umap.plot)
dev.off()
```
![scRNA_atlas](https://github.com/ucdavis/quonlab/blob/master/development/deconvolution/deconvAllen/deconvTutorial/mouse_scRNA_atlas.png)

## Deconv. model setup
Now that we have our cell type labels and marker gene set finalized it's time to build the deconvolution model.
`<NAME>` requires 3 inputs:
```text
scRNA Atlas                   cell x gene
scRNA Atlas cell type labels  cell x 1
mixture                       sample x gene
```

Furthermore, we have to determine which genes to initially pass into our model. In this case we are going to use the common
variable and marker genes across mouse and mouse.patchseq.

```R
## The common set of variable and marker genes which our cell type models will reconstruct
deconv.gene.set = Reduce(union, list(union(marker.genes, VariableFeatures(mouse)), VariableFeatures(mouse.patchseq)))
deconv.gene.set = intersect(deconv.gene.set, rownames(mouse))
deconv.gene.set = intersect(deconv.gene.set, rownames(mouse.patchseq))

## Lets define a mask for the marker genes (1=marker, 0=non_marker)
marker_gene_mask = rep(0, length(deconv.gene.set))
marker_gene_mask[which(deconv.gene.set %in% marker.genes)] = 1

## Gather the scRNA atlas and subset to genes
mouse.data    = as.matrix(GetAssayData(mouse, "scale.data")[deconv.gene.set,])

## Gather the PatchSeq gene expr. and subset to genes
patchseq.data = as.matrix(GetAssayData(mouse.patchseq, "scale.data"))[deconv.gene.set,]

## Now we create the deconvolution model!
deconvModel = deconv$deconvModel(component_data  = t(mouse.data),                             ## Altas
                                 component_label = as.array(as.character(mouse@meta.data$celltype.deconv)),   ## Cell type annotations for each cell in Atlas
                                 mixture_data    = t(patchseq.data))                          ## Mixtures
```

## Deconvolution of mouse PatchSeq gene expr.
Now, we take our defined `deconvModel` and run the deconvolution task!
```R
## Now lets run deconvolution, I have left the hyperparameters exposed here so you can see the different ways that the model could be tuned.
deconvModel$deconvolve(## Masks
                      marker_gene_mask = marker_gene_mask,
                      ## Training steps
                      max_steps_component  = 1500L,
                      max_steps_proportion = 25000L,
                      max_steps_combined   = 5000L,
                      ## Learning rate
                      component_learning_rate  = 1e-3,
                      proportion_learning_rate = 1e-3,
                      combined_learning_rate   = 1e-3,
                      ## Early stopping
                      early_stopping       = 'False',
                      early_min_step       = 500L,
                      max_patience         = 100L,
                      ## Output
                      log_step             = 25L,
                      print_step           = 100L,
                      ## Seed
                      seed                 = 1234L,
                      ## VAE params
                      KL_weight            = 1.0,
                      num_latent_dims      = 32L,
                      num_layers           = 3L,
                      hidden_unit_2power   = 9L,
                      decoder_var          = "per_sample",
                      ## Batch size
                      batch_size_component = 150L,
                      batch_size_mixture   = 300L,
                       ## Training method
                       training_method     = 'train', ## train: use all scRNA data, "valid": split scRNA into training and testing sets
                       ## Deconvolution flags
                       batch_correction    = 'False', ## Don't change this for now.
                       decay_lr            = 'False',
                       batch_norm_layers   = 'True',
                       corr_check          = 'False',
                       log_results         = 'True',
                       log_samples         = 'False',
                       cuda_device         = 0L)
```

## Getting the deconv results
While `deconvModel` holds all the input data and results, it is a python class which R can't interpret. So we have a helper function to convert `deconvModel` to an R S4 class which will contain all the results from deconvolution plus data annotations.
```R
## Convert the python class to an R S4 class. check: str(deconvResults)
deconvResults = convertDeconv(deconvModel,                 ## Deconv model object
                              deconv.gene.set,             ## Features used during deconvolution
                              colnames(mouse.data),        ## scRNA cell ids
                              colnames(mouse.patchseq))    ## mixture sample ids

## Get the final proportion estimates
proportions = deconvResults@proportions$final

## Call cell type from proportions
est.labels = deconvModel$celltypes[apply(proportions, 1, which.max)]

## Plot the estimated proportions
pdf(paste0("./mouse_patchseq_deconv_summary.pdf"), width=14, height=12)
propHeatmap(proportions, mixture.labels=est.labels)
# model_corrCheck(deconvResults, deconvModel)         ## Check model
# model_tsneViz(deconvResults, deconvModel)           ## Plot results
dev.off()

## Save the deconvolution model and the marker gene set
save(deconvResults, marker.anno.df, file=paste0("./deconv_tutorial_results.rda"))
```
## Proportion sanity check
Let's take a quick look at the proportions and make sure the marker gene expression supports the estimated cell type proportions.
```R
proportionSanityCheck(deconvResults,
                      mouse.patchseq,
                      data = "scale.data", ## "data" or "scale.data" from Seurat object
                      filename="./mouse_patchseq_deconv_sanity.png",
                      marker.anno=marker.anno.df)
```
![proportion_check](https://github.com/ucdavis/quonlab/blob/master/development/deconvolution/deconvAllen/deconvTutorial/mouse_patchseq_example_run.png)

## Purification sanity check (IN PROGRESS)
<!-- Finally, `<NAME>` also outputs all the mixture samples purified to each cell type in the scRNAseq atlas. Let's take a look at the purified and scRNA data.
```R
## Let's first check all the types of purified data we have!
print(names(deconvResults@deconv_data$purified$train))

## We extract the purified data for the `Sst` cell type as follows:
pure.data.sst = t(deconvResults@deconv_data$purified$train$Sst) ## Make gene x cell for now
rownames(pure.data.sst) = deconvResults@genes
colnames(pure.data.sst) = deconvResults@mixtureCellNames

## Now lets plot the scRNA and purified data together to see how our model performed.
## First step, convert the purified data to a Seurat object (in a somewhat hacky way)
pure.seurat = CreateSeuratObject(pure.data.sst)
pure.seurat@assays$RNA@data = pure.seurat@assays$RNA@scale.data = pure.data.sst

## Merge the pure and scRNA Data
merged.seurat = merge(mouse, pure.seurat)
merged.seurat@assays$RNA@scale.data = cbind(GetAssayData(mouse, "scale.data")[rownames(pure.seurat),],
                                            GetAssayData(pure.seurat, "scale.data")[rownames(pure.seurat),])

## Update the celltype.deconv for the Sst purified data
merged.seurat@meta.data$celltype.deconv = c(mouse@meta.data$celltype.deconv, rep("Sst", ncol(pure.seurat)))
merged.seurat@meta.data$datatype = c(rep("scRNA", ncol(mouse)), rep("purifiedSst", ncol(pure.seurat)))

## Run dimensionality reduction on just the marker genes
merged.seurat <- RunPCA(merged.seurat, npcs = 30, features=marker.genes, verbose = FALSE)
merged.seurat <- RunUMAP(merged.seurat, reduction = "pca", dims = 1:30)

## Plot
celltype.plot <- DimPlot(merged.seurat, reduction = "umap", group.by = "celltype.deconv", label=TRUE, repel=TRUE)
datatype.plot <- DimPlot(merged.seurat, reduction = "umap", group.by = "datatype", label=TRUE, repel=TRUE)

## Save our plot
library(cowplot)
png("./mouse_pure_scRNA_atlas.png", width=24, height=12, res=150,  units="in")
plot_grid(celltype.plot, datatype.plot)
dev.off()
```
![purified_check]() -->

```R
sessionInfo()


```
