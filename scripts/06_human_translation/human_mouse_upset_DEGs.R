##############################
## Upset plot of KEGG terms ##
##############################
# Author : Osman Sharifi

#############################
## Load required libraries ##
#############################
library(UpSetR)
library(dplyr)
library(enrichR)
library(openxlsx)
library(glue)
library(ggplot2)
library(purrr)
library(rrvgo)
library(data.table)
library(stringr)
library(org.Hs.eg.db)
library(forcats)
library(Hmisc)
library(ggsci)

#load data
human_degs <- read.csv('/Users/osman/Documents/GitHub/PCB_mouse_snRNAseq_BACKUP_20251031/09_human_DEG_analysis/limmaVoomCC/Human_rett_PCBvs_noPCB_sig_DEGs.csv')

human_glut_list <- human_degs %>% filter(Cell_type == "Glutamatergic")
human_gaba_list <- human_degs %>% filter(Cell_type == "GABAergic")
human_nonn_list <- human_degs %>% filter(Cell_type == "Non-Neuronal")
mouse_degs <- read.csv('/Users/osman/Documents/GitHub/PCB_mouse_snRNAseq_BACKUP_20251031/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE.csv')
# Convert the SYMBOL column to uppercase
mouse_degs <- mouse_degs %>%
  mutate(SYMBOL = toupper(SYMBOL))
mouse_glut_list <- mouse_degs %>% filter(Cell_type == "Glutamatergic")
mouse_gaba_list <- mouse_degs %>% filter(Cell_type == "GABAergic")
mouse_nonn_list <- mouse_degs %>% filter(Cell_type == "Non-neuronal")

#######################
## Create upset plot ##
#######################
listInput <- list(Glutamatergic_human = human_glut_list$SYMBOL, 
                  Glutamatergic_mouse = mouse_glut_list$SYMBOL,
                  GABAergic_human = human_gaba_list$SYMBOL, 
                  GABAergic_mouse = mouse_gaba_list$SYMBOL,
                  Non_neuronal_human = human_nonn_list$SYMBOL,
                  Non_neuronal_mouse = mouse_nonn_list$SYMBOL)
# Find common SYMBOLs across all lists
common_symbols <- Reduce(intersect, listInput)

# View the result
common_symbols
# Find common SYMBOLs between GABAergic_human and GABAergic_mouse
common_gaba_symbols <- intersect(listInput$GABAergic_human, listInput$GABAergic_mouse)
common_gaba_symbols <- stringr::str_remove(common_gaba_symbols, "\\d$")
newlist <- list(GABAergic_human = human_nonn_list$SYMBOL,  GABAergic_mouse = mouse_nonn_list$SYMBOL)
# View the result
common_gaba_symbols <- Reduce(intersect, newlist)
common_gaba_symbols

base_path <- "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/09_human_DEG_analysis/limmaVoomCC/"
pdf(glue("{base_path}upset_human_mouse_DEGs.pdf"))
upset(fromList(listInput), sets = c('Glutamatergic_human','Glutamatergic_mouse', 'GABAergic_human', 'GABAergic_mouse', 'Non_neuronal_human', 'Non_neuronal_mouse'), keep.order = TRUE)
dev.off()

#############################
## Define slimGO function  ##
#############################

slimGO <- function(GO = GO,
                   tool = c("enrichR", "rGREAT", "GOfuncR"),
                   annoDb = annoDb,
                   plots = FALSE,
                   threshold = 0.7){
  
  if(tool == "enrichR"){
    GO <- GO %>%
      data.table::rbindlist(idcol = "Gene Ontology") %>%
      dplyr::as_tibble() %>%
      dplyr::filter(grepl("^GO_", `Gene Ontology`)) %>%
      dplyr::mutate(Term = stringr::str_extract(.$Term, "\\(GO.*")) %>%
      dplyr::mutate(Term = stringr::str_replace_all(.$Term, "[//(//)]", "")) %>%
      dplyr::mutate("Gene Ontology" = dplyr::case_when(
        grepl("Biological_Process", `Gene Ontology`) ~ "BP",
        grepl("Cellular_Component", `Gene Ontology`) ~ "CC",
        grepl("Molecular_Function", `Gene Ontology`) ~ "MF")) %>%
      dplyr::select(p = P.value, go = Term, "Gene Ontology") %>%
      dplyr::filter(p < 0.05) %>%
      dplyr::filter(!is.na(go) & go != "")
  }
  
  print(glue::glue("Submitting results from {tool} to rrvgo..."))
  print(glue::glue("GO terms per ontology:"))
  print(table(GO$`Gene Ontology`))
  
  .slim <- function(GO, ont, annoDb, plots, tool, threshold){
    
    GO_filtered <- GO %>% dplyr::filter(`Gene Ontology` == ont)
    
    if(nrow(GO_filtered) < 2){
      print(glue::glue("Skipping {ont}: fewer than 2 significant terms"))
      return(NULL)
    }
    
    print(glue::glue("rrvgo is now slimming {ont} GO terms from {tool} ({nrow(GO_filtered)} terms)"))
    
    tryCatch({
      simMatrix <- rrvgo::calculateSimMatrix(GO_filtered$go, orgdb = annoDb, ont = ont, method = "Rel")
      
      if(is.null(simMatrix) || nrow(simMatrix) < 2){
        print(glue::glue("Skipping {ont}: similarity matrix too small"))
        return(NULL)
      }
      
      reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, setNames(-log10(GO_filtered$p), GO_filtered$go), threshold = threshold, orgdb = annoDb)
      
      if(plots == TRUE){ 
        p <- rrvgo::scatterPlot(simMatrix, reducedTerms)
        plot(p)
        rrvgo::treemapPlot(reducedTerms) 
      }
      
      print(glue::glue("There are {max(reducedTerms$cluster)} clusters in your GO {ont} terms from {tool}"))
      return(dplyr::as_tibble(reducedTerms))
      
    }, error = function(e){
      print(glue::glue("Error in {ont}: {e$message}"))
      return(NULL)
    })
  }
  
  ontologies <- GO %>%
    dplyr::select(`Gene Ontology`) %>%
    table() %>%
    names()
  
  slimmed_list <- ontologies %>%
    purrr::set_names() %>%
    purrr::map(~.slim(GO = GO, ont = ., annoDb = annoDb, tool = tool, plots = plots, threshold = threshold))
  
  slimmed_list <- slimmed_list[!sapply(slimmed_list, is.null)]
  
  if(length(slimmed_list) == 0){
    stop("No ontologies could be processed successfully")
  }
  
  slimmed <- slimmed_list %>%
    dplyr::bind_rows(.id = "Gene Ontology") %>%
    dplyr::inner_join(GO, by = c("Gene Ontology", "go")) %>%
    dplyr::filter(term == as.character(parentTerm)) %>%
    dplyr::mutate("-log10.p-value" = -log10(p)) %>%
    dplyr::mutate("Gene Ontology" = dplyr::recode_factor(`Gene Ontology`, 
                                                         "BP" = "Biological Process", 
                                                         "CC" = "Cellular Component", 
                                                         "MF" = "Molecular Function")) %>%
    dplyr::arrange(dplyr::desc(`-log10.p-value`)) %>%
    dplyr::select("Gene Ontology", Term = term, "-log10.p-value")
  
  return(slimmed)
}

#############################
## Define GOplot function  ##
#############################

GOplot <- function(slimmedGO){
  slimmedGO %>%
    dplyr::group_by(`Gene Ontology`) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Term = stringr::str_trim(Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
    dplyr::mutate(Term = stringr::str_wrap(Term, 50)) %>%
    dplyr::mutate(Term = factor(Term, levels = unique(Term[order(forcats::fct_rev(`Gene Ontology`), `-log10.p-value`)]))) %>%
    ggplot2::ggplot(ggplot2::aes(x = Term, y = `-log10.p-value`, fill = `Gene Ontology`, group = `Gene Ontology`)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "Black") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](p))) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14), 
                   axis.title.y = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 14), 
                   legend.title = ggplot2::element_text(size = 14),
                   strip.text = ggplot2::element_text(size = 14))
}

#############################
## Load and process DEGs   ##
#############################

# Load mouse DEGs, convert to uppercase, and remove trailing digit (only if not after a decimal)
mouse_degs <- read.csv('/Users/osman/Documents/GitHub/PCB_mouse_snRNAseq_BACKUP_20251031/07_mosiacism/2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE.csv') %>%
  mutate(SYMBOL = toupper(SYMBOL)) %>%
  mutate(SYMBOL = stringr::str_remove(SYMBOL, "(?<!\\.)\\d$"))

# Load human DEGs, remove trailing digit (only if not after a decimal), and standardize cell type names
human_degs <- read.csv('/Users/osman/Documents/GitHub/PCB_mouse_snRNAseq_BACKUP_20251031/09_human_DEG_analysis/limmaVoomCC/Human_rett_PCBvs_noPCB_sig_DEGs.csv') %>%
  mutate(SYMBOL = stringr::str_remove(SYMBOL, "(?<!\\.)\\d$")) %>%
  mutate(Cell_type = case_when(
    Cell_type == "Non-Neuronal" ~ "Non-neuronal",
    TRUE ~ Cell_type
  ))

# Verify the fix
print("Sample human symbols after fix:")
print(head(human_degs$SYMBOL, 20))
print("Sample mouse symbols after fix:")
print(head(mouse_degs$SYMBOL, 20))

# Check for any symbols ending with . (should be none)
print("Human symbols ending with '.':")
print(human_degs$SYMBOL[grepl("\\.$", human_degs$SYMBOL)] %>% head(10))
print("Mouse symbols ending with '.':")
print(mouse_degs$SYMBOL[grepl("\\.$", mouse_degs$SYMBOL)] %>% head(10))

# Check cell types
print("Human cell types:")
print(table(human_degs$Cell_type))
print("Mouse cell types:")
print(table(mouse_degs$Cell_type))

# Create cell type lists
human_glut_list <- human_degs %>% filter(Cell_type == "Glutamatergic")
human_gaba_list <- human_degs %>% filter(Cell_type == "GABAergic")
human_nonn_list <- human_degs %>% filter(Cell_type == "Non-neuronal")

mouse_glut_list <- mouse_degs %>% filter(Cell_type == "Glutamatergic")
mouse_gaba_list <- mouse_degs %>% filter(Cell_type == "GABAergic")
mouse_nonn_list <- mouse_degs %>% filter(Cell_type == "Non-neuronal")

# Find common genes between human and mouse
common_gaba_symbols <- intersect(human_gaba_list$SYMBOL, mouse_gaba_list$SYMBOL)
common_glut_symbols <- intersect(human_glut_list$SYMBOL, mouse_glut_list$SYMBOL)
common_nonn_symbols <- intersect(human_nonn_list$SYMBOL, mouse_nonn_list$SYMBOL)

# Verify intersections
print(glue("Human GABAergic DEGs: {nrow(human_gaba_list)}"))
print(glue("Mouse GABAergic DEGs: {nrow(mouse_gaba_list)}"))
print(glue("Common GABAergic genes: {length(common_gaba_symbols)}"))

print(glue("Human Glutamatergic DEGs: {nrow(human_glut_list)}"))
print(glue("Mouse Glutamatergic DEGs: {nrow(mouse_glut_list)}"))
print(glue("Common Glutamatergic genes: {length(common_glut_symbols)}"))

print(glue("Human Non-neuronal DEGs: {nrow(human_nonn_list)}"))
print(glue("Mouse Non-neuronal DEGs: {nrow(mouse_nonn_list)}"))
print(glue("Common Non-neuronal genes: {length(common_nonn_symbols)}"))

print("Sample common GABAergic genes:")
print(head(common_gaba_symbols, 10))

#############################
## Run enrichR analysis    ##
#############################

directory_path <- "/Users/osman/Documents/GitHub/PCB_mouse_snRNAseq_BACKUP_20251031/enrichR_results"
if (!dir.exists(directory_path)) dir.create(directory_path, recursive = TRUE)

# Define databases
dbs <- c("GO_Biological_Process_2023",
         "GO_Molecular_Function_2023",
         "GO_Cellular_Component_2023",
         "KEGG_2021_Human",
         "Panther_2016",
         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")

#############################
## GABAergic enrichment    ##
#############################

if(length(common_gaba_symbols) >= 5){
  tryCatch({
    print("Running GABAergic enrichment...")
    
    enriched_results <- enrichR::enrichr(common_gaba_symbols, dbs)
    
    print("Results per database:")
    print(sapply(enriched_results, nrow))
    
    # Save raw results
    enriched_results %>%
      purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %>%
      openxlsx::write.xlsx(file = glue::glue("{directory_path}/GABAergic_enrichr.xlsx"))
    
    # Run slimGO and plot
    slimmed_results <- enriched_results %>%
      purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %>%
      slimGO(tool = "enrichR", annoDb = "org.Hs.eg.db", plots = FALSE)
    
    openxlsx::write.xlsx(slimmed_results, file = glue::glue("{directory_path}/GABAergic_rrvgo_enrichr.xlsx"))
    
    go_plot <- GOplot(slimmed_results)
    
    ggplot2::ggsave(glue::glue("{directory_path}/GABAergic_enrichr_plot.pdf"),
                    plot = go_plot,
                    height = 8.5,
                    width = 10)
    
    print("GABAergic enrichment complete!")
    
  }, error = function(e) {
    print(glue::glue("ERROR in GABAergic enrichment: {e$message}"))
  })
} else {
  print(glue("Only {length(common_gaba_symbols)} common GABAergic genes found - skipping enrichment (need at least 5)"))
}

#############################
## Glutamatergic enrichment ##
#############################

if(length(common_glut_symbols) >= 5){
  tryCatch({
    print("Running Glutamatergic enrichment...")
    
    enriched_results <- enrichR::enrichr(common_glut_symbols, dbs)
    
    print("Results per database:")
    print(sapply(enriched_results, nrow))
    
    enriched_results %>%
      purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %>%
      openxlsx::write.xlsx(file = glue::glue("{directory_path}/Glutamatergic_enrichr.xlsx"))
    
    slimmed_results <- enriched_results %>%
      purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %>%
      slimGO(tool = "enrichR", annoDb = "org.Hs.eg.db", plots = FALSE)
    
    openxlsx::write.xlsx(slimmed_results, file = glue::glue("{directory_path}/Glutamatergic_rrvgo_enrichr.xlsx"))
    
    go_plot <- GOplot(slimmed_results)
    
    ggplot2::ggsave(glue::glue("{directory_path}/Glutamatergic_enrichr_plot.pdf"),
                    plot = go_plot,
                    height = 8.5,
                    width = 10)
    
    print("Glutamatergic enrichment complete!")
    
  }, error = function(e) {
    print(glue::glue("ERROR in Glutamatergic enrichment: {e$message}"))
  })
} else {
  print(glue("Only {length(common_glut_symbols)} common Glutamatergic genes found - skipping enrichment (need at least 5)"))
}

#############################
## Non-neuronal enrichment ##
#############################

if(length(common_nonn_symbols) >= 5){
  tryCatch({
    print("Running Non-neuronal enrichment...")
    
    enriched_results <- enrichR::enrichr(common_nonn_symbols, dbs)
    
    print("Results per database:")
    print(sapply(enriched_results, nrow))
    
    enriched_results %>%
      purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %>%
      openxlsx::write.xlsx(file = glue::glue("{directory_path}/Non-neuronal_enrichr.xlsx"))
    
    slimmed_results <- enriched_results %>%
      purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %>%
      slimGO(tool = "enrichR", annoDb = "org.Hs.eg.db", plots = FALSE)
    
    openxlsx::write.xlsx(slimmed_results, file = glue::glue("{directory_path}/Non-neuronal_rrvgo_enrichr.xlsx"))
    
    go_plot <- GOplot(slimmed_results)
    
    ggplot2::ggsave(glue::glue("{directory_path}/Non-neuronal_enrichr_plot.pdf"),
                    plot = go_plot,
                    height = 8.5,
                    width = 10)
    
    print("Non-neuronal enrichment complete!")
    
  }, error = function(e) {
    print(glue::glue("ERROR in Non-neuronal enrichment: {e$message}"))
  })
} else {
  print(glue("Only {length(common_nonn_symbols)} common Non-neuronal genes found - skipping enrichment (need at least 5)"))
}

print("All enrichment analyses complete!")
