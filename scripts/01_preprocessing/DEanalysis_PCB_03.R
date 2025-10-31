### Differential Expression Analysis of scRNA-seq Data - 00 - Splitting Seurat Object ###

# Splitting Seurat Object by Activated and Unactivated Neurons

# set up 

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

s.obj.name = "all_female_P30" #change this to the name of the Seurat object you're working with

load(glue::glue("/Users/karineier/Documents/Mecp2/scRNA-seq/{s.obj.name}.RData"))

s.obj = all_female_P30 #change this to the name of the Seurat object you're working with

dir.create(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/{s.obj.name}"))

setwd(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/{s.obj.name}"))

### Split Seurat object into activated and unactivated neurons ###

activated.neurons.df = ifelse((s.obj$predicted.id == "L2_3_IT" | 
                                 s.obj$predicted.id == "L4" | 
                                 s.obj$predicted.id == "L5" | 
                                 s.obj$predicted.id == "L6" |
                                 s.obj$predicted.id ==  "Pvalb" | 
                                 s.obj$predicted.id == "Lamp5" | 
                                 s.obj$predicted.id == "Vip" | 
                                 s.obj$predicted.id == "Sst" | 
                                 s.obj$predicted.id == "Sncg") &
                                s.obj@assays$RNA@counts["Fos",]>0, "-activated", "")

activated.names = sapply(1:length(s.obj$predicted.id), function(x){
  new.id = gsub("$", glue::glue("{activated.neurons.df[[x]]}"), s.obj$predicted.id[[x]])
  return(new.id)
})

s.obj$cell.type = activated.names

Idents(s.obj) = "cell.type"

save(s.obj, activated.neurons.df, activated.names, file="DEanalysis_00.RData")
