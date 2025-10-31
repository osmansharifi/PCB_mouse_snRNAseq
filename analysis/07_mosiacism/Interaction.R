###############
## libraries ##
###############
library(ggplot2)
library(Seurat)
library(dplyr)

###############
## load data ##
###############
load("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/PEBBLES_parsed.RData")
xist_data_total <- FetchData(PEBBLES_soupx, vars = c("Xist", "Group", "Mecp2_allele", "broad_class", "Genotype", "Treatment"))

#############################
## Create interaction plot ##
#############################
xist_data_total <- xist_data_total %>%
  filter(Mecp2_allele != "Unparsed")

xist_data <- xist_data_total %>%
  filter(Genotype != "WT")
         
xist_data <- xist_data %>%
  mutate(Mecp2_allele = factor(Mecp2_allele, levels = c("WT_Mecp2", "MUT_Mecp2")))

# Perform three-way ANOVA (including broad_class)
anova_result <- aov(Xist ~ Mecp2_allele * Treatment * broad_class, data = xist_data)
anova_summary <- summary(anova_result)

# Extract F-values and p-values
anova_stats <- data.frame(
  Term = rownames(anova_summary[[1]]),
  F_value = anova_summary[[1]]$`F value`,
  P_value = anova_summary[[1]]$`Pr(>F)`
)

# Format ANOVA results for annotation
anova_text <- paste(
  "Three-Way ANOVA:\n",
  "Mecp2_allele: F = ", round(anova_stats$F_value[1], 2), ", p = ", round(anova_stats$P_value[1], 4), "\n",
  "Treatment: F = ", round(anova_stats$F_value[2], 2), ", p = ", round(anova_stats$P_value[2], 4), "\n",
  "broad_class: F = ", round(anova_stats$F_value[3], 2), ", p = ", round(anova_stats$P_value[3], 4), "\n",
  "Mecp2_allele:Treatment: F = ", round(anova_stats$F_value[4], 2), ", p = ", round(anova_stats$P_value[4], 4), "\n",
  "Mecp2_allele:broad_class: F = ", round(anova_stats$F_value[5], 2), ", p = ", round(anova_stats$P_value[5], 4), "\n",
  "Treatment:broad_class: F = ", round(anova_stats$F_value[6], 2), ", p = ", round(anova_stats$P_value[6], 4), "\n",
  "Mecp2_allele:Treatment:broad_class: F = ", round(anova_stats$F_value[7], 2), ", p = ", round(anova_stats$P_value[7], 4),
  sep = ""
)
# Create the interaction plot
ggplot(xist_data, aes(x = Mecp2_allele, y = Xist, color = broad_class, shape = Treatment, group = interaction(Treatment, broad_class))) +
           stat_summary(fun = mean, geom = "line", linewidth = 1.5) +  # Use `linewidth` for lines
           stat_summary(fun = mean, geom = "point", size = 5) +      # Add points for means
           theme_minimal() +
           labs(title = "Interaction Plot: Xist Expression by Mecp2 Allele and Treatment",
                x = "Mecp2_allele",
                y = "Xist Expression Level",
                color = "Cell Type",
                shape = "Treatment") +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 90, size = 18, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 18, face = 'bold', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    # Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"), # Text size
    title = element_text(size = 18, face = "bold")
  )
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/interaction_plot.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

# Create the interaction plot with ANOVA statistics
ggplot(xist_data, aes(x = Mecp2_allele, y = Xist, color = broad_class, shape = Treatment, group = interaction(Treatment, broad_class))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +  # Use `linewidth` for lines
  stat_summary(fun = mean, geom = "point", size = 5) +      # Add points for means
  theme_minimal() +
  labs(title = "Interaction Plot: Xist Expression by Mecp2 Allele and Treatment",
       x = "Mecp2_allele",
       y = "Xist Expression Level",
       color = "Cell Type",
       shape = "Treatment") +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 90, size = 18, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 18, face = 'bold', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 18, face = "bold"),
    title = element_text(size = 18, face = "bold")
  ) +
  annotate("text", x = Inf, y = Inf, label = anova_text, hjust = 1.1, vjust = 1.1, size = 5, color = "black")


### Testing

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(patchwork)  # For arranging plots side by side

# Filter data for MUT and WT samples
xist_data_total <- xist_data_total %>%
  filter(Mecp2_allele != "Unparsed")

xist_data_mut <- xist_data_total %>%
  filter(Genotype != "WT")

xist_data_wt <- xist_data_total %>%
  filter(Genotype == "WT")

# Ensure Mecp2_allele is a factor with the correct levels
xist_data_mut <- xist_data_mut %>%
  mutate(Mecp2_allele = factor(Mecp2_allele, levels = c("WT_Mecp2", "MUT_Mecp2")))

# Perform three-way ANOVA for MUT samples (as in your original code)
anova_result <- aov(Xist ~ Mecp2_allele * Treatment * broad_class, data = xist_data_mut)
anova_summary <- summary(anova_result)

# Extract F-values and p-values
anova_stats <- data.frame(
  Term = rownames(anova_summary[[1]]),
  F_value = anova_summary[[1]]$`F value`,
  P_value = anova_summary[[1]]$`Pr(>F)`
)

# Format ANOVA results for annotation
anova_text <- paste(
  "Three-Way ANOVA:\n",
  "Mecp2_allele: F = ", round(anova_stats$F_value[1], 2), ", p = ", round(anova_stats$P_value[1], 4), "\n",
  "Treatment: F = ", round(anova_stats$F_value[2], 2), ", p = ", round(anova_stats$P_value[2], 4), "\n",
  "broad_class: F = ", round(anova_stats$F_value[3], 2), ", p = ", round(anova_stats$P_value[3], 4), "\n",
  "Mecp2_allele:Treatment: F = ", round(anova_stats$F_value[4], 2), ", p = ", round(anova_stats$P_value[4], 4), "\n",
  "Mecp2_allele:broad_class: F = ", round(anova_stats$F_value[5], 2), ", p = ", round(anova_stats$P_value[5], 4), "\n",
  "Treatment:broad_class: F = ", round(anova_stats$F_value[6], 2), ", p = ", round(anova_stats$P_value[6], 4), "\n",
  "Mecp2_allele:Treatment:broad_class: F = ", round(anova_stats$F_value[7], 2), ", p = ", round(anova_stats$P_value[7], 4),
  sep = ""
)

# Create the interaction plot for MUT samples
plot_mut <- ggplot(xist_data_mut, aes(x = Mecp2_allele, y = Xist, color = broad_class, shape = Treatment, group = interaction(Treatment, broad_class))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  theme_minimal() +
  labs(title = "MUT Samples: Xist Expression by Mecp2 Allele and Treatment",
       x = "Mecp2_allele",
       y = "Xist Expression Level",
       color = "Cell Type",
       shape = "Treatment") +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    plot.title = element_text(size = 18, face = 'bold', vjust = 1),
    axis.text.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.text.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    legend.text = element_text(size = 18, face = "bold"),
    title = element_text(size = 18, face = "bold")
  )

# Create the interaction plot for WT samples
plot_wt <- ggplot(xist_data_wt, aes(x = Mecp2_allele, y = Xist, color = broad_class, shape = Treatment, group = interaction(Treatment, broad_class))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  theme_minimal() +
  labs(title = "WT Samples: Xist Expression by Mecp2 Allele and Treatment",
       x = "Mecp2_allele",
       y = "Xist Expression Level",
       color = "Cell Type",
       shape = "Treatment") +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    plot.title = element_text(size = 18, face = 'bold', vjust = 1),
    axis.text.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.text.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    legend.text = element_text(size = 18, face = "bold"),
    title = element_text(size = 18, face = "bold")
  )

# Arrange plots side by side using patchwork
combined_plot <- plot_mut + plot_wt

# Display the combined plot
combined_plot
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/interaction_plot_MUT_WT.pdf",
                device = NULL,
                height = 8.5,
                width = 12)


## Testing

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(patchwork)  # For arranging plots side by side

# Filter data for MUT and WT samples
xist_data_total <- xist_data_total %>%
  filter(Mecp2_allele != "Unparsed")

xist_data_mut <- xist_data_total %>%
  filter(Genotype != "WT")

xist_data_wt <- xist_data_total %>%
  filter(Genotype == "WT")

# Ensure Mecp2_allele is a factor with the correct levels
xist_data_mut <- xist_data_mut %>%
  mutate(Mecp2_allele = factor(Mecp2_allele, levels = c("WT_Mecp2", "MUT_Mecp2")))

# Perform three-way ANOVA for MUT samples (as in your original code)
anova_result <- aov(Xist ~ Mecp2_allele * Treatment * broad_class, data = xist_data_mut)
anova_summary <- summary(anova_result)

# Extract F-values and p-values
anova_stats <- data.frame(
  Term = rownames(anova_summary[[1]]),
  F_value = anova_summary[[1]]$`F value`,
  P_value = anova_summary[[1]]$`Pr(>F)`
)

# Format ANOVA results for annotation
anova_text <- paste(
  "Three-Way ANOVA:\n",
  "Mecp2_allele: F = ", round(anova_stats$F_value[1], 2), ", p = ", round(anova_stats$P_value[1], 4), "\n",
  "Treatment: F = ", round(anova_stats$F_value[2], 2), ", p = ", round(anova_stats$P_value[2], 4), "\n",
  "broad_class: F = ", round(anova_stats$F_value[3], 2), ", p = ", round(anova_stats$P_value[3], 4), "\n",
  "Mecp2_allele:Treatment: F = ", round(anova_stats$F_value[4], 2), ", p = ", round(anova_stats$P_value[4], 4), "\n",
  "Mecp2_allele:broad_class: F = ", round(anova_stats$F_value[5], 2), ", p = ", round(anova_stats$P_value[5], 4), "\n",
  "Treatment:broad_class: F = ", round(anova_stats$F_value[6], 2), ", p = ", round(anova_stats$P_value[6], 4), "\n",
  "Mecp2_allele:Treatment:broad_class: F = ", round(anova_stats$F_value[7], 2), ", p = ", round(anova_stats$P_value[7], 4),
  sep = ""
)

# Create the interaction plot for MUT samples
plot_mut <- ggplot(xist_data_mut, aes(x = Mecp2_allele, y = Xist, color = broad_class, shape = Treatment, group = interaction(Treatment, broad_class))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  theme_minimal() +
  labs(title = "MUT Samples: Xist Expression by Mecp2 Allele and Treatment",
       x = "Mecp2_allele",
       y = "Xist Expression Level",
       color = "Cell Type",
       shape = "Treatment") +
  scale_y_continuous(limits = c(3.4, ceiling(max(xist_data_total$Xist, na.rm = TRUE))),  # Y-axis starts at 3.4
                     breaks = seq(3.4, ceiling(max(xist_data_total$Xist, na.rm = TRUE)))) +  # Increments by 0.1
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    plot.title = element_text(size = 18, face = 'bold', vjust = 1),
    axis.text.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.text.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    legend.text = element_text(size = 18, face = "bold"),
    title = element_text(size = 18, face = "bold")
  )

# Create the interaction plot for WT samples
plot_wt <- ggplot(xist_data_wt, aes(x = Mecp2_allele, y = Xist, color = broad_class, shape = Treatment, group = interaction(Treatment, broad_class))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  theme_minimal() +
  labs(title = "WT Samples: Xist Expression by Mecp2 Allele and Treatment",
       x = "Mecp2_allele",
       y = "Xist Expression Level",
       color = "Cell Type",
       shape = "Treatment") +
  scale_y_continuous(limits = c(3.4, ceiling(max(xist_data_total$Xist, na.rm = TRUE))),  # Y-axis starts at 3.4
                     breaks = seq(3.4, ceiling(max(xist_data_total$Xist, na.rm = TRUE)))) +  # Increments by 0.1
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    plot.title = element_text(size = 18, face = 'bold', vjust = 1),
    axis.text.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.text.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    legend.text = element_text(size = 18, face = "bold"),
    title = element_text(size = 18, face = "bold")
  )

# Arrange plots side by side using patchwork
combined_plot <- plot_mut + plot_wt

# Display the combined plot
combined_plot

# Boxplot for Xist by Mecp2_allele and Treatment
boxplot(Xist ~ Mecp2_allele * Treatment, data = xist_data_mut,
        col = c("lightblue", "lightgreen"), main = "Xist by Mecp2_allele and Treatment",
        xlab = "Mecp2_allele:Treatment", ylab = "Xist")

# Boxplot for Xist by Mecp2_allele and broad_class
boxplot(Xist ~ Mecp2_allele * broad_class, data = xist_data_mut,
        col = c("lightpink", "lightyellow"), main = "Xist by Mecp2_allele and Broad Class",
        xlab = "Mecp2_allele:Broad Class", ylab = "Xist")

# interaction plot of genotype and treatment
xist_data <- xist_data_total %>%
  mutate(Genotype = factor(Genotype),  # Ensure Genotype is a factor
         Treatment = factor(Treatment))  # Ensure Treatment is a factor

#############################
## Perform Two-Way ANOVA   ##
#############################
# Two-way ANOVA for Xist ~ Genotype * Treatment
anova_result <- aov(Xist ~ Genotype * Treatment, data = xist_data_total)
anova_summary <- summary(anova_result)

# Extract F-values and p-values
anova_stats <- data.frame(
  Term = rownames(anova_summary[[1]]),
  F_value = anova_summary[[1]]$`F value`,
  P_value = anova_summary[[1]]$`Pr(>F)`
)

# Format ANOVA results for annotation
anova_text <- paste(
  "Two-Way ANOVA:\n",
  "Genotype: F = ", round(anova_stats$F_value[1], 2), ", p = ", round(anova_stats$P_value[1], 4), "\n",
  "Treatment: F = ", round(anova_stats$F_value[2], 2), ", p = ", round(anova_stats$P_value[2], 4), "\n",
  "Genotype:Treatment: F = ", round(anova_stats$F_value[3], 2), ", p = ", round(anova_stats$P_value[3], 4),
  sep = ""
)

#############################
## Create Interaction Plot ##
#############################
ggplot(xist_data_total, aes(x = Genotype, y = Xist, color = Treatment, group = Treatment)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +  # Lines for means
  stat_summary(fun = mean, geom = "point", size = 5) +       # Points for means
  theme_minimal() +
  labs(title = "Interaction Plot: Xist Expression by Genotype and Treatment",
       x = "Genotype",
       y = "Xist Expression Level",
       color = "Treatment") +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 90, size = 18, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 18, face = 'bold', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    # Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"), # Text size
    title = element_text(size = 18, face = "bold")
  ) +
  annotate("text", x = Inf, y = Inf, label = anova_text, hjust = 1.1, vjust = 1.1, size = 5, color = "black")  # Add ANOVA results as annotation

# Save the plot
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/interaction_plot_genotype_treatment.pdf",
                plot = interaction_plot,
                device = NULL,
                height = 8.5,
                width = 12)
