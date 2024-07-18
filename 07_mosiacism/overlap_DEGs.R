######################
## Overlap analysis ##
######################
library(magrittr)
library(VennDiagram)
library(grDevices)
library(dplyr)
library(Seurat)
library(glue)
library(scCustomize)
library(edgeR)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(GeneOverlap)
library(enrichR)
library(readr)
library(stringr)
library(UpSetR)
################
## load files ##
################
# Set the working directory to the parent folder containing the subdirectories
setwd("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/07_mosiacism/")

# List of file paths
file_paths <- c(
  "1_AllcellsVsAllcells_from_MUTPCB_WTPCB/sig_DEGs_1_AllcellsVsAllcells_from_MUTPCB_WTPCB.csv",
  "2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_2_AllcellsVsAllcells_from_MUTPCB_MUTVEHICLE.csv",
  "3_WTcellsVsWTcells_from_WTPCB_WTVEHICLE/sig_DEGs_3_WTcellsVsWTcells_from_WTPCB_WTVEHICLE.csv",
  "4_WTcellsVsWTcells_from_MUTPCB_WTPCB/sig_DEGs_4_WTcellsVsWTcells_from_MUTPCB_WTPCB.csv",
  "5_WTcellsVsWTcells_from_MUTPCB_WTVEHICLE/sig_DEGs_5_WTcellsVsWTcells_from_MUTPCB_WTVEHICLE.csv",
  "6_AllcellsVsAllcells_from_MUTVEHICLE_WTVEHICLE/sig_DEGs_6_AllcellsVsAllcells_from_MUTVEHICLE_WTVEHICLE.csv",
  "7_WTcellsVsWTcells_from_MUTVEHICLE_WTVEHICLE/sig_DEGs_7_WTcellsVsWTcells_from_MUTVEHICLE_WTVEHICLE.csv",
  "8_WTcellsVsWTcells_from_MUTVEHICLE_WTPCB/sig_DEGs_8_WTcellsVsWTcells_from_MUTVEHICLE_WTPCB.csv",
  "9_WTcellsVsWTcells_from_HETPCB_HETVEHICLE/sig_DEGs_9_WTcellsVsWTcells_from_HETPCB_HETVEHICLE.csv",
  "10_MUTcellsVsMUTcells_from_HETPCB_HETVEHICLE/sig_DEGs_10_MUTcellsVsMUTcells_from_HETPCB_HETVEHICLE.csv",
  "11_WTcellsVsMUTcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_11_WTcellsVsMUTcells_from_MUTPCB_MUTVEHICLE.csv",
  "12_MUTcellsVsWTcells_from_MUTPCB_MUTVEHICLE/sig_DEGs_12_MUTcellsVsWTcells_from_MUTPCB_MUTVEHICLE.csv",
  "13_MUTcellsVsWTcells_within_MUTVEHICLE/sig_DEGs_13_MUTcellsVsWTcells_within_MUTVEHICLE.csv",
  "14_MUTcellsVsWTcells_within_MUTPCB/sig_DEGs_14_MUTcellsVsWTcells_within_MUTPCB.csv",
  "15_MUTcellsVsWTcells_from_MUTPCB_WTPCB/sig_DEGs_15_MUTcellsVsWTcells_from_MUTPCB_WTPCB.csv"
)

# Read each file into a list of data frames
exp_list <- lapply(file_paths, read.csv, header = TRUE)

# Extract the filenames without extension
file_names <- tools::file_path_sans_ext(basename(file_paths))

# Add DEG_experiment column to each data frame
for (i in seq_along(exp_list)) {
  exp_list[[i]]$DEG_experiment <- file_names[i]
}

combined_df <- do.call(rbind, exp_list)

# Write the concatenated data frame to a new CSV file
write_csv(combined_df, "Summary_data/all_experiments_sig_DEGs.csv")

###############################
## Number of DEGs up or down ##
###############################
summary_table <- combined_df %>%
  group_by(DEG_experiment, Cell_type) %>%
  summarise(
    UP = sum(logFC > 0, na.rm = TRUE),
    DOWN = sum(logFC < 0, na.rm = TRUE)
  ) %>%
  ungroup()

print(summary_table)
write_csv(summary_table, "Summary_data/all_experiments_sig_DEGs.csv")


# Ensure summary_table is in the right format
summary_table <- summary_table %>%
  gather(key = "logFC_type", value = "count", UP, DOWN)

summary_table$Regulation = summary_table$logFC_type

summary_table <- summary_table %>%
  mutate(
    # Extract the number between underscores and convert it to numeric
    OrderNumber = as.numeric(str_extract(DEG_experiment, "(?<=_)[0-9]+(?=_)"))
  ) %>%
  arrange(OrderNumber) %>% # Arrange the rows based on the extracted number
  select(-OrderNumber) # Optionally remove the OrderNumber column
summary_table$DEG_experiment <- factor(summary_table$DEG_experiment, 
                                       levels = rev(unique(summary_table$DEG_experiment)))

summary_table <- summary_table %>%
  filter(Cell_type != "Non-neuronal")

ggplot(summary_table, aes(x = DEG_experiment, y = count, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() + # Flip the coordinates to make the bars horizontal
  facet_wrap(~ Cell_type) +
  theme_minimal() +
  labs(x = "DEG Experiment", y = "DEG count", fill = "Regulation") +
  scale_fill_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  scale_y_continuous(labels = abs) + # Use absolute values for labels
  theme(
    text = element_text(size=16),
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 16, face = 'bold'),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")
  )
ggplot2::ggsave(("Summary_data/all_experiments_logFC_DEGs_non-neuronalremoved.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

#################
## UPset plots ##
#################

aggregate_symbols <- function(df) {
  df %>%
    group_by(DEG_experiment) %>%
    summarise(SYMBOL = list(unique(SYMBOL))) %>%
    pull(SYMBOL)
}

glutamatergic_sets <- aggregate_symbols(combined_df %>%
                                          filter(Cell_type == "Glutamatergic"))

gabaergic_sets <- aggregate_symbols(combined_df %>%
                                      filter(Cell_type == "GABAergic"))

combined_df %>%
  ggplot(aes(x = SYMBOL)) +
  geom_bar() +
  facet_wrap(~ DEG_experiment, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "SYMBOL", y = "Count", title = "Overlap of SYMBOL for each DEG_experiment")
