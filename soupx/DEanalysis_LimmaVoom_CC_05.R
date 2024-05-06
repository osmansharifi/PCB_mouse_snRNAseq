###################
## Load Packages ##
###################
packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Glimma", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis", "RColorBrewer", "BiocParallel", "DEsingle", "enrichR", "DMRichR", "org.Mm.eg.db", "AnnotationDbi")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

enrichR:::.onAttach()

############################
## Load Pathway Functions ##
############################
# SlimGO
slimGO <- function(GO = GO,
                   tool = c("enrichR", "rGREAT", "GOfuncR"),
                   annoDb = annoDb,
                   plots = FALSE,
                   threshold = 0.7){
  
  if(tool == "enrichR"){
    GO <- GO %>%
      data.table::rbindlist(idcol = "Gene Ontology") %>%
      dplyr::as_tibble() %>%
      dplyr::filter(`Gene Ontology` %in% c("GO_Biological_Process_2018",
                                           "GO_Cellular_Component_2018",
                                           "GO_Molecular_Function_2018")) %>% 
      dplyr::mutate(Term = stringr::str_extract(.$Term, "\\(GO.*")) %>%
      dplyr::mutate(Term = stringr::str_replace_all(.$Term, "[//(//)]",""), "") %>%
      dplyr::mutate("Gene Ontology" = dplyr::case_when(`Gene Ontology` == "GO_Biological_Process_2018" ~ "BP",
                                                       `Gene Ontology` == "GO_Cellular_Component_2018" ~ "CC",
                                                       `Gene Ontology` == "GO_Molecular_Function_2018" ~ "MF")) %>%
      dplyr::select(p = P.value, go = Term, "Gene Ontology") %>% 
      dplyr::filter(p < 0.05)
    
  }else if(tool == "rGREAT"){
    
    GO <-  GO %>%
      data.table::rbindlist(idcol = "Gene Ontology") %>%
      dplyr::as_tibble() %>%
      dplyr::mutate("Gene Ontology" = dplyr::case_when(`Gene Ontology` == "GO Biological Process" ~ "BP",
                                                       `Gene Ontology` == "GO Cellular Component" ~ "CC",
                                                       `Gene Ontology` == "GO Molecular Function" ~ "MF")) %>%
      dplyr::select(p = Hyper_Raw_PValue, go = ID, "Gene Ontology") %>% 
      dplyr::filter(p < 0.05)
    
  }else if(tool == "GOfuncR"){
    
    GO <- GO$results %>%
      dplyr::as_tibble() %>%
      dplyr::mutate("Gene Ontology" = dplyr::case_when(ontology == "biological_process" ~ "BP",
                                                       ontology == "cellular_component" ~ "CC",
                                                       ontology == "molecular_function" ~ "MF")) %>%
      dplyr::select(p = raw_p_overrep, go = node_id, "Gene Ontology") %>% 
      dplyr::filter(p < 0.05)
    
  }else{
    stop(glue("{tool} is not supported, please choose either enrichR, rGREAT, or GOfuncR [Case Sensitive]"))
  }
  
  print(glue::glue("Submiting results from {tool} to rrvgo..."))
  
  .slim <- function(GO = GO,
                    ont = ont,
                    annoDb = annoDb,
                    plots = plots,
                    tool = tool,
                    threshold = threshold){
    GO <- GO %>%
      dplyr::filter(`Gene Ontology` == ont)
    
    print(glue::glue("rrvgo is now slimming {ont} GO terms from {tool}"))
    
    simMatrix  <- rrvgo::calculateSimMatrix(GO$go,
                                            orgdb = annoDb,
                                            ont = ont,
                                            method = "Rel")
    
    reducedTerms <- rrvgo::reduceSimMatrix(simMatrix,
                                           setNames(-log10(GO$p), GO$go),
                                           threshold = threshold,
                                           orgdb = annoDb) 
    
    if(plots == TRUE){
      p <- rrvgo::scatterPlot(simMatrix, reducedTerms) # Doesn't plot otherwise
      plot(p) 
      rrvgo::treemapPlot(reducedTerms)
    }
    
    print(glue::glue("There are {max(reducedTerms$cluster)} clusters in your GO {ont} terms from {tool}"))
    
    reducedTerms %>%   
      dplyr::as_tibble() %>%
      return()
  }
  
  slimmed <- GO %>%
    dplyr::select(`Gene Ontology`) %>%
    table() %>%
    names() %>% 
    purrr::set_names() %>%
    purrr::map_dfr(~.slim(GO = GO,
                          ont = .,
                          annoDb = annoDb,
                          tool = tool,
                          plots = plots,
                          threshold = threshold),
                   .id = "Gene Ontology") %>%
    dplyr::inner_join(GO) %>%
    dplyr::filter(term == as.character(parentTerm)) %>%
    dplyr::mutate("-log10.p-value" = -log10(p)) %>%
    dplyr::mutate("Gene Ontology" = dplyr::recode_factor(`Gene Ontology`,
                                                         "BP" = "Biological Process",
                                                         "CC" = "Cellular Component",
                                                         "MF" = "Molecular Function")) %>%
    dplyr::arrange(dplyr::desc(`-log10.p-value`)) %>% 
    dplyr::select("Gene Ontology", Term = term, "-log10.p-value") %>% 
    return()
}

# GOplot
GOplot <- function(slimmedGO = slimmedGO){
  
  slimmedGO %>% 
    dplyr::group_by(`Gene Ontology`) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_wrap(.$Term, 50)) %>% 
    #dplyr::mutate(Term = stringr::str_trunc(.$Term, 45, side = "right")) %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$`Gene Ontology`), .$`-log10.p-value`)]))) %>% 
    ggplot2::ggplot(ggplot2::aes(x = Term,
                                 y = `-log10.p-value`,
                                 fill = `Gene Ontology`,
                                 group = `Gene Ontology`)) +
    ggplot2::geom_bar(stat = "identity",
                      position = ggplot2::position_dodge(),
                      color = "Black") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](p))) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                   axis.title.y = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 14),
                   legend.title = ggplot2::element_text(size = 14),
                   strip.text = ggplot2::element_text(size = 14)) %>% 
    return()
}

#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

###############
## Load Data ##
###############
#setwd(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rettcort_labeled"))
load("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rettcort_labeled/DEanalysis_01.1.RData")
s.obj = human #change this to the name of the Seurat object you started with

dir.create("limmaVoomCC")
setwd("limmaVoomCC")

##############################
## DEG and Pathway Analysis ##
##############################

#DGEList_high = DGEList_high[-c(which(names(DGEList_high)!="Micro"))] #removing non-neuronal cells due to lack of sufficient cell numbers across genotypes

cell_types_all = names(DGEList_high)

#cell_types_all = cell_types_all[-(c(which(grepl("Sncg-activated", cell_types_all)=="TRUE")))] #removing Sncg-activated from analysis bc not enough cells

cell_types_all = cell_types_all[c(1:15)] # not enough cells for these cell types

for (i in cell_types_all) {
  DGEList_high[[i]]$samples$group = ifelse(grepl("WT", colnames(DGEList_high[[i]]$counts))=="TRUE", "WT", "RTT")
}

cell_types_all = cell_types_all[c(1:15)]

for (i in cell_types_all) {
  
  dir.create(glue::glue("{i}"))
  dir.create(glue::glue("{i}/QC"))
  dir.create(glue::glue("{i}/plotData"))
  
  print(glue::glue("Normalizing {i} cells"))
  
  #subsetting design matrix by cell type
  design$genotype = factor(design$genotype, levels=c("CTRL", "RTT"))
  
  # Raw density of log-CPM values
  
  L <- mean(DGEList_high[[i]]$samples$lib.size)*1e-6
  M <- median(DGEList_high[[i]]$samples$lib.size)*1e-6
  
  logCPM <- cpm(DGEList_high[[i]], log=TRUE)
  nsamples <- ncol(DGEList_high[[i]])
  col <- polychrome_palette
  
  pdf(glue::glue("{i}/density_plot.pdf"), height=8.5, width=5.5)
  
  plot(density(logCPM[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="Density Plot of Highly Expressed Genes", xlab="Log-cpm")
  for(j in 2:nsamples){
    den <- density(logCPM[,j])
    lines(den$x, den$y, col=col[j], lwd=2)
  }
  legend("topright", colnames(DGEList_high[[i]]$counts), text.col=col, bty="n", cex=0.5)
  
  dev.off()
  
  #Glimma::glMDSPlot(DGEList_high[[i]],
                    #groups = design$genotype,
                   # path = getwd(),
                   # folder = "interactivePlots",
                   # html = glue::glue("{i}_MDS-Plot"),
                   # launch = FALSE)
  
  
  mm <- model.matrix(~genotype + cell_cycle + percent.mito,
                     data = design)
  
  
  # Voom
  pdf(glue::glue("{i}/QC/voom_mean-variance_trend.pdf"), height=8.5, width=11)
  voomLogCPM <- voom(DGEList_high[[i]],
                     mm,
                     plot=TRUE)
  dev.off()
  
  correlations <- duplicateCorrelation(voomLogCPM,
                                       mm,
                                       block = design$sample_ID)
  
  correlations <- correlations$consensus.correlation
  
  # Boxplots of logCPM values before and after voom normalization
  pdf(glue::glue("{i}/QC/normalization_boxplots.pdf"), height=8.5, width=11)
  par(mfrow=c(1,2))
  
  boxplot(logCPM, las=2, col=col, main="")
  title(main="A. Unnormalized data", ylab="Log-cpm")
  
  boxplot(voomLogCPM$E, las=2, col=col, main="")
  title(main="B. Normalized data", ylab="Log-cpm")
  
  dev.off()
  
  # Fitting linear models in limma
  
  print(glue::glue("Testing {i} cells for differential expression"))
  
  fit <- lmFit(voomLogCPM,
               mm,
               correlation = correlations,
               block = design$sample_ID)
  
  head(coef(fit))
  
  # Save normalized expression values
  voomLogCPM$E %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "geneSymbol") %>%
    openxlsx::write.xlsx(glue::glue("{i}/plotData/voomLogCPM.xlsx"))
  
  # Create DEG List 
  
  print(glue::glue("Creating DEG list of {i} cells for genotype"))
  
  efit <- fit %>%
    contrasts.fit(coef="genotypeMUTANT") %>%
    eBayes()
  
  pdf(glue::glue("{i}/QC/final_model_mean-variance_trend.pdf"), height=8.5, width=11)
  
  plotSA(efit, main=glue::glue("Final model {i}: Mean-variance trend"))
  
  dev.off()
  
  Glimma::glimmaMA(efit,
                   dge = DGEList_high[[i]],
                   path = getwd(),
                   html = glue::glue("interactivePlots/{i}_MDA-Plot.html"),
                   launch = FALSE)
  
  # Top DEGs
  
  DEGs <- efit %>%
    topTable(sort.by = "P", n = Inf) %>%
    tibble::rownames_to_column() %>%
    tibble::as_tibble() %>%
    dplyr::rename(SYMBOL = rowname) %>%
    dplyr::mutate(FC = dplyr::case_when(logFC >0 ~ 2^logFC,
                                        logFC <0 ~ -1/(2^logFC))) %>%
    dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %T>%
    openxlsx::write.xlsx(file=glue::glue("{i}/DEGs.xlsx")) %>%
    dplyr::filter(P.Value < 0.05) %T>%
    openxlsx::write.xlsx(file=glue::glue("{i}/sig_DEGs.xlsx"))
  
  print(glue::glue("GO and Pathway analysis of {i} cells"))
  
  enrichR:::.onAttach()
  
  tryCatch({
    DEGs %>% 
      dplyr::select(SYMBOL) %>%
      purrr::flatten() %>%
      enrichR::enrichr(c("GO_Biological_Process_2018",
                         "GO_Molecular_Function_2018",
                         "GO_Cellular_Component_2018",
                         "KEGG_2019_Mouse",
                         "Panther_2016",
                         "Reactome_2016",
                         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
      purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
      openxlsx::write.xlsx(file=glue::glue("{i}/enrichr.xlsx")) %>%
      slimGO(tool = "enrichR",
                      annoDb = "org.Mm.eg.db",
                      plots = FALSE) %T>%
      openxlsx::write.xlsx(file = glue::glue("{i}/rrvgo_enrichr.xlsx")) %>%
      GOplot() %>%
      ggplot2::ggsave(glue::glue("{i}/enrichr_plot.pdf"),
                      plot = .,
                      device = NULL,
                      height = 8.5,
                      width = 10) },
    error = function(error_condition) {
      print(glue::glue("ERROR: Gene Ontology pipe did not finish for {i} cells"))
    })
  print(glue::glue("The pipeline has finished for {i} cells"))
  
}
