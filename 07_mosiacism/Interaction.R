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
           stat_summary(fun = mean, geom = "line", linewidth = 1) +  # Use `linewidth` for lines
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


# Create the interaction plot with ANOVA statistics
ggplot(xist_data, aes(x = Mecp2_allele, y = Xist, color = broad_class, shape = Treatment, group = interaction(Treatment, broad_class))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +  # Use `linewidth` for lines
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
