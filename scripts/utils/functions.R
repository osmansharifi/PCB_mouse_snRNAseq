#' ============================================================================
#' Utility Functions for PCB snRNA-seq Analysis
#' ============================================================================

#' Load and Validate Seurat Object
load_seurat <- function(file_path, required_metadata = NULL) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  seurat_obj <- readRDS(file_path)
  
  if (!is(seurat_obj, "Seurat")) {
    stop("Loaded object is not a Seurat object")
  }
  
  if (!is.null(required_metadata)) {
    missing_cols <- setdiff(required_metadata, colnames(seurat_obj@meta.data))
    if (length(missing_cols) > 0) {
      stop(paste("Missing metadata columns:", paste(missing_cols, collapse = ", ")))
    }
  }
  
  return(seurat_obj)
}

#' Calculate QC Metrics
calculate_qc_metrics <- function(seurat_obj, 
                                 mt_pattern = "^MT-",
                                 rb_pattern = "^RB") {
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = rb_pattern)
  
  return(seurat_obj)
}

#' Apply QC Filters
apply_qc_filters <- function(seurat_obj,
                             nCount_min = 500,
                             nCount_max = 50000,
                             nFeature_min = 200,
                             percent_mt_max = 5,
                             verbose = TRUE) {
  
  n_before <- ncol(seurat_obj)
  
  seurat_obj <- subset(
    seurat_obj,
    subset = (nCount_RNA >= nCount_min &
              nCount_RNA <= nCount_max &
              nFeature_RNA >= nFeature_min &
              percent.mt <= percent_mt_max)
  )
  
  n_after <- ncol(seurat_obj)
  
  if (verbose) {
    cat(sprintf("QC filtering: %d → %d cells (%.1f%% retained)\n",
                n_before, n_after, n_after/n_before*100))
  }
  
  return(seurat_obj)
}

#' Normalize and Scale Data
normalize_and_scale <- function(seurat_obj,
                               normalization_method = "LogNormalize",
                               scale_factor = 10000,
                               verbose = TRUE) {
  
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = verbose
  )
  
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = verbose
  )
  
  seurat_obj <- ScaleData(seurat_obj, verbose = verbose)
  
  return(seurat_obj)
}

#' Run DEG Analysis for All Cell Types
run_deg_all_celltypes <- function(seurat_obj,
                                 cell_type_col = "cell_type",
                                 ident1,
                                 ident2,
                                 group_by = "treatment",
                                 min_log2fc = 0.25,
                                 min_pct = 0.1) {
  
  cell_types <- unique(seurat_obj@meta.data[[cell_type_col]])
  deg_list <- list()
  
  for (cell_type in cell_types) {
    subset_obj <- subset(seurat_obj, subset = !!sym(cell_type_col) == cell_type)
    
    deg_results <- tryCatch({
      FindMarkers(
        subset_obj,
        ident.1 = ident1,
        ident.2 = ident2,
        group.by = group_by,
        test.use = "wilcox",
        min.log2fc = min_log2fc,
        min.pct = min_pct,
        verbose = FALSE
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(deg_results) && nrow(deg_results) > 0) {
      deg_results$gene <- rownames(deg_results)
      deg_results$cell_type <- cell_type
      deg_list[[cell_type]] <- deg_results
    }
  }
  
  deg_all <- bind_rows(deg_list)
  return(deg_all)
}

#' Save Seurat Object
save_seurat <- function(seurat_obj, file_path) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(seurat_obj, file = file_path)
}

#' Check if Files Exist
check_files_exist <- function(file_paths) {
  missing <- !file.exists(file_paths)
  
  if (any(missing)) {
    missing_files <- file_paths[missing]
    stop(paste("Required files not found:\n", 
               paste("  -", missing_files, collapse = "\n")))
  }
  
  invisible(TRUE)
}
