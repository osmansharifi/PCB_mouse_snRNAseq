##### Summary of scRNA-seq DEG Analysis #####

library(openxlsx)
library(glue)
library(tidyr)
library(ggplot2)
library(viridus)

cell_types = c("Astro", "L2_3_IT", "L4", "L5", "L6", "Lamp5", "Non-neuronal", "Oligo", "Pvalb", "Sncg", "Sst", "Vip")

list <- vector(mode="list", length=length(cell_types))

list <- lapply(cell_types, function(cellType) {
  openxlsx::read.xlsx(glue::glue("/Users/karineier/Documents/scRNA-seq/Genotype/noDreamWeights/Females/{cellType}/DEGs.xlsx"), rowNames = TRUE)
})

names(list) = cell_types

numDEGs = sapply(cell_types, function(cellType) {
  length(which(list[[cellType]][,4] < 0.1))
})

plotData = data.frame(cell_type = cell_types, numDEGs = as.numeric(numDEGs))

plot = ggplot(plotData, aes(x=cell_type, y=numDEGs, fill=numDEGs)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(legend.position="NULL", axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_viridis(option="plasma") +
  xlab("") +
  ylab("") 

ggsave("Summary_Genotype_DEGs_by_celltype_females_noDreamWeights.pdf", width=8.5, height=11)




