#' ============================================================================
#' Project: PCB Exposure and Gene Expression in Mouse Brain (snRNA-seq)
#' Script: SoupX Ambient RNA Correction
#' Purpose: Remove ambient RNA contamination using SoupX
#' 
#' Author: Osman Sharifi
#' Date Created: 2024-10-31
#' Last Modified: 2024-10-31
#' 
#' Input:  Cellranger output from 10X encapsulation
#' Output: data/processed/seurat_soupx_corrected.rds
#'         Seurat object with corrected counts
#' 
#' Dependencies:
#'   - Seurat (v4.3.0+)
#'   - SoupX
#'   - functions from: scripts/utils/functions.R
#'   - logging from: scripts/utils/logging.R
#' 
#' Notes:
#'   - SoupX contamination fraction auto-estimated from data
#'   - Uses relative paths for portability
#'   - Random seed set to 42 for reproducibility
#' ============================================================================

# ============ 1. SETUP & PARAMETERS ============
library(Seurat)
library(SoupX)
library(dplyr)

# Source utility functions
source("scripts/utils/functions.R")
source("scripts/utils/logging.R")

set.seed(42)

# Define parameters
SOUPX_PARAMS <- list(
  min_cells = 10,
  min_features = 200,
  description = "SoupX ambient RNA correction parameters"
)

log_message("Starting SoupX ambient RNA correction", level = "INFO")

# ============ 2. DEFINE SAMPLE INFORMATION ============

# List of samples to process
samples <- c(
  "24_PCB_WT",
  "25_VEHICLE_WT",
  "27_PCB_HET_2",
  "27_PCB_HET",
  "28_VEHICLE_HET",
  "29_VEHICLE_WT",
  "30_VEHICLE_WT_2",
  "30_VEHICLE_WT",
  "31_PCB_WT",
  "36_PCB_HET",
  "37_PCB_WT_2",
  "37_PCB_WT",
  "38_VEHICLE_HET",
  "39_PCB_HET",
  "40_VEHICLE_HET_2",
  "40_VEHICLE_HET"
)

# Update this path to your actual cellranger output directory
# Using relative path from repository root
cellranger_path <- "data/raw/cellranger_output/"

log_message(
  sprintf("Processing %d samples from: %s", length(samples), cellranger_path),
  level = "INFO"
)

# ============ 3. HELPER FUNCTION: PROCESS SAMPLE WITH SOUPX ============

process_sample_soupx <- function(sample_name, folder_path) {
  
  tryCatch({
    # Load 10X data
    sample_path <- paste0(folder_path, sample_name, "/outs/")
    
    if (!dir.exists(sample_path)) {
      log_message(
        sprintf("WARNING: Path not found for %s: %s", sample_name, sample_path),
        level = "WARNING"
      )
      return(NULL)
    }
    
    # Load and estimate contamination
    sc <- load10X(sample_path)
    sc <- autoEstCont(sc)
    
    # Get contamination fraction
    contam_fraction <- sc$metaData$nUMIs[2] / sc$metaData$nUMIs[1]
    log_message(
      sprintf("%s: Contamination fraction = %.3f", sample_name, contam_fraction),
      level = "INFO"
    )
    
    # Adjust counts
    corrected_counts <- adjustCounts(sc)
    
    # Add sample name prefix to barcodes
    colnames(corrected_counts) <- paste0(sample_name, "_", colnames(corrected_counts))
    
    return(corrected_counts)
    
  }, error = function(e) {
    log_message(
      sprintf("ERROR processing %s: %s", sample_name, e$message),
      level = "ERROR"
    )
    return(NULL)
  })
}

# ============ 4. PROCESS ALL SAMPLES ============

log_message("Processing samples with SoupX...", level = "INFO")

# Process each sample
corrected_matrices <- list()

for (i in seq_along(samples)) {
  log_progress(i, length(samples), prefix = "Processing:")
  
  corrected <- process_sample_soupx(samples[i], cellranger_path)
  
  if (!is.null(corrected)) {
    corrected_matrices[[samples[i]]] <- corrected
  }
}

log_message(
  sprintf("Successfully processed %d/%d samples", 
          length(corrected_matrices), length(samples)),
  level = "SUCCESS"
)

# ============ 5. COMBINE AND CREATE SEURAT OBJECT ============

log_message("Combining samples and creating Seurat object...", level = "INFO")

# Combine all count matrices
combined_counts <- do.call(cbind, corrected_matrices)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = combined_counts,
  project = "PEBBLES_SoupX",
  min.cells = SOUPX_PARAMS$min_cells,
  min.features = SOUPX_PARAMS$min_features,
  names.field = 2,
  names.delim = "_"
)

log_message(
  sprintf("Created Seurat object: %d cells, %d genes",
          ncol(seurat_obj), nrow(seurat_obj)),
  level = "SUCCESS"
)

# ============ 6. ADD METADATA ============

log_message("Adding sample metadata...", level = "INFO")

# Extract sample names from barcodes
sample_names <- sub("_[^_]+$", "", colnames(seurat_obj))

# Add metadata columns
seurat_obj$sample <- sample_names
seurat_obj$genotype <- ifelse(grepl("WT", sample_names), "WT",
                              ifelse(grepl("HET", sample_names), "HET", NA))
seurat_obj$treatment <- ifelse(grepl("PCB", sample_names), "PCB",
                               ifelse(grepl("VEHICLE", sample_names), "VEHICLE", NA))

log_message("Metadata added successfully", level = "SUCCESS")

# ============ 7. SAVE RESULTS ============

log_message("Saving corrected Seurat object...", level = "INFO")

# Create output directory
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# Save as RDS
output_file <- "data/processed/seurat_soupx_corrected.rds"
saveRDS(seurat_obj, file = output_file)

log_file_save(
  output_file,
  description = "SoupX-corrected Seurat object"
)

# ============ 8. SESSION INFO ============

log_session_info()

log_message("Script complete!", level = "SUCCESS")
