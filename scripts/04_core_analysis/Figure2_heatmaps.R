library(pheatmap)
library(tidyr)
library(dplyr)
library(gridExtra)
DEGs <- read.csv("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_DEGs_filtered.csv")

## top DEGs
top_genes <- DEGs %>%
  group_by(DEG_test, Directory) %>%
  arrange(desc(logFC)) %>%
  slice(1:2) %>%  # Top 10 upregulated genes
  bind_rows(
    DEGs %>%
      group_by(DEG_test, Directory) %>%
      arrange(logFC) %>%
      slice(1:2)  # Top 10 downregulated genes
  ) %>%
  ungroup() %>%
  distinct(SYMBOL, DEG_test, Directory, .keep_all = TRUE)

# Spread the data into a matrix format
heatmap_data <- top_genes %>%
  select(SYMBOL, logFC, DEG_test, Directory) %>%
  pivot_wider(
    names_from = c(DEG_test, Directory), 
    values_from = logFC, 
    values_fill = list(logFC = 0)
  )

# Convert to a matrix
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$SYMBOL

# Generate column metadata
split_names <- strsplit(colnames(heatmap_matrix), "_")
col_metaData <- data.frame(
  Comparison = sapply(split_names, `[`, 1),
  Genotype = sapply(split_names, `[`, 3),
  CellType = sapply(split_names, `[`, 4)
)
rownames(col_metaData) <- colnames(heatmap_matrix)

# Define color schemes
comparison_colors <- c("WTvsHET" = "#EE7E30", "VehicleVsPCB" = "#5D9AD3")
genotype_colors <- c("Vehicle" = "#98D352", "WTs" = "#FF7F0E", "PCBs" = "#1F77B4", "HETs" = "#D62728")
#celltype_colors <- c("Non-neuronal" = "#8D80BA", "L4" = "#DD902A", "Lamp5" = "#ED9DC7", "L2" = "#EA5B57", "L5" = "#F4D642", "Sncg" = "#1541BC", "Sst" = "#0098C1", "Vip" = "#5CEABE", "Pvalb" = "#2F9349", "L6" = "#2ED852", "Oligo" = "#C77CFF")
celltype_colors <- c("L2" = "#EA5B57", "L4" = "#DD902A","L5" = "#F4D642", "L6" = "#2ED852", "Sst" = "#0098C1", "Vip" = "#5CEABE", "Pvalb" = "#2F9349", "Sncg" = "#1541BC", "Oligo" = "#C77CFF","Lamp5" = "#ED9DC7", "Non-neuronal" = "#8D80BA")
col <- list(
  Comparison = comparison_colors,
  Genotype = genotype_colors,
  CellType = celltype_colors
)

# Define comparisons
comparisons <- unique(col_metaData$Comparison)

# Generate the heatmaps
heatmap_list <- lapply(comparisons, function(comp) {
  subset_matrix <- heatmap_matrix[, col_metaData$Comparison == comp]
  subset_col_metaData <- col_metaData[col_metaData$Comparison == comp, ]
  
  # Remove rows with zero variance
  row_variances <- apply(subset_matrix, 1, var, na.rm = TRUE)
  subset_matrix <- subset_matrix[row_variances > 0, ]
  
  # Remove columns with zero variance
  col_variances <- apply(subset_matrix, 2, var, na.rm = TRUE)
  subset_matrix <- subset_matrix[, col_variances > 0]
  
  # Debugging: Print the dimensions and summary of subset_matrix
  cat("Comparison:", comp, "\n")
  cat("Dimensions of subset_matrix:", dim(subset_matrix), "\n")
  print(summary(subset_matrix))
  
  pheatmap(
    subset_matrix,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = subset_col_metaData,
    annotation_colors = col,
    show_rownames = TRUE,  # Show gene names on the left
    show_colnames = FALSE,
    color = colorRampPalette(c("#7A9DEA", "white", "#F44780"))(100),
    legend_labels = list(color = "logFC"),  # Label the color bar
    silent = TRUE
  )
})

# Save the heatmaps side by side
pdf("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/testing_heatmap.pdf", width = 16, height = 8)

# Use grid.arrange with adjusted layout
grid.arrange(
  grobs = lapply(heatmap_list, function(x) x[[4]]), 
  ncol = 2, 
  widths = c(1, 1)  # Equal widths for both heatmaps
)

dev.off()

