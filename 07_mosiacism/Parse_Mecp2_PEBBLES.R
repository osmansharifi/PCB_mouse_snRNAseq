# This program adds Mecp2_allele counts back to the PEBBLES seurat object

####################
## load libraries ##
####################
packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "dwtools", "devtools", "dplyr", "patchwork", "scCustomize")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

########################
## load Seurat object ##
########################
setwd("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/")
load("../06_hdWGCNA/PEBBLES_clean.RData")

#################################
## Prepare Mecp2 allele counts ##
#################################
counts_table <- read.table("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/PEBBLES_Mecp2_alles.txt", sep="\t", header=FALSE)
#headers for the txt file
names(counts_table) <- c("Barcode", "UMI", "WT", "MUT", "Body", "Sample")
head(counts_table)
#counts_table <- counts_table[grepl("WT", counts_table$Sample), ]
counts_table <- counts_table[!duplicated(counts_table$Barcode), ]
counts_table <- counts_table %>%
  filter(!(grepl("WT", Sample) & MUT > 0))
Mecp2_wt_mut_counts = counts_table

##############################################
## Integrate counts into the PEBBLES object ##
##############################################
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
rownames(PEBBLES_soupx@assays$RNA@counts) = c(rownames(PEBBLES_soupx@assays$RNA@counts)[-c(18611:18612)], "Mecp2-WT", "Mecp2-MUT")

# checking to make sure they are added
PEBBLES_soupx@assays$RNA@counts[18611:18612, 1:5]
PEBBLES_soupx@meta.data <- PEBBLES_soupx@meta.data %>%
  mutate(Mecp2_allele = case_when(
    WT_Mecp2 > 0 ~ "WT_Mecp2",
    MUT_Mecp2 > 0 ~ "MUT_Mecp2",
    TRUE ~ "Unparsed"
  ))

##########################
## Plot the parsed data ##
##########################
FeaturePlot_scCustom(seurat_object = PEBBLES_soupx, features = c('WT_Mecp2','MUT_Mecp2'), split.by = "Genotype",slot = 'counts', pt.size = 0.7)
ggplot2::ggsave("Mecp2_parsed_clean.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
Idents(PEBBLES_soupx) <- "cell_type"
DimPlot_scCustom(seurat_object = PEBBLES_soupx, split.by = "Mecp2_allele", pt.size = 0.7, order = TRUE)


