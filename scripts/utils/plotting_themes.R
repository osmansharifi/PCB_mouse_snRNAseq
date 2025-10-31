#' ============================================================================
#' ggplot2 Themes and Color Palettes for PCB snRNA-seq Analysis
#' ============================================================================

library(ggplot2)

#' PCB Analysis ggplot2 Theme
theme_pcb <- function() {
  theme_minimal() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      color = "gray40",
      margin = margin(b = 10)
    ),
    
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "gray30", size = 0.3),
    
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "gray50"),
    
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5)
  )
}

#' Cell Type Color Palette
cell_type_palette <- function() {
  list(
    "Excitatory" = "#E74C3C",
    "Inhibitory" = "#3498DB",
    "Astrocyte" = "#F39C12",
    "Oligodendrocyte" = "#9B59B6",
    "OPC" = "#1ABC9C",
    "Microglia" = "#E67E22",
    "Endothelial" = "#16A085",
    "Pericyte" = "#27AE60"
  )
}

#' Treatment Color Palette
treatment_palette <- function() {
  list(
    "PCB" = "#E74C3C",
    "Vehicle" = "#95A5A6",
    "Control" = "#95A5A6"
  )
}

#' Direction Color Palette (Up/Down)
direction_palette <- function() {
  list(
    "up" = "#E74C3C",
    "down" = "#3498DB"
  )
}

#' Create volcano plot
plot_volcano <- function(deg_results,
                        fc_threshold = 0.25,
                        p_threshold = 0.05,
                        title = "Volcano Plot",
                        label_top_n = 5) {
  
  deg_results <- deg_results %>%
    mutate(
      direction = case_when(
        avg_log2FC > fc_threshold & p_val_adj < p_threshold ~ "up",
        avg_log2FC < -fc_threshold & p_val_adj < p_threshold ~ "down",
        TRUE ~ "ns"
      )
    )
  
  p <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = direction), size = 2.5, alpha = 0.7) +
    scale_color_manual(
      values = c(
        "up" = "#E74C3C",
        "down" = "#3498DB",
        "ns" = "#BDC3C7"
      ),
      name = "Direction"
    ) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "gray50", alpha = 0.5) +
    theme_pcb() +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value"
    )
  
  return(p)
}
