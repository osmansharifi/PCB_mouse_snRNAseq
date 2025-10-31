library(UpSetR)
library(dplyr)
library(ggplot2)

# Get a list of all directories
all_directories <- list.dirs(path = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx", full.names = TRUE)

# Create an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each directory
for (dir in all_directories) {
  # Get the file path for enrichr.xlsx in the current directory
  keggs_file <- file.path(dir, "enrichr.xlsx")
  
  # Check if the file exists
  if (file.exists(keggs_file)) {
    # Read the enrichr.xlsx file
    keggs_data <- read_excel(keggs_file, sheet = "KEGG_2019_Mouse")
    
    # Add a new column with the directory name
    keggs_data$Directory <- basename(dir)
    
    # Get the parent directory name
    parent_dir <- dirname(dir)
    
    # Add a new column with the parent directory name
    keggs_data$DEG_test <- basename(parent_dir)
    
    # Combine with existing data
    combined_data <- rbind(combined_data, keggs_data)
  }
}

# unfiltered total keggs
write.csv(combined_data, file = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_keggs_unfiltered.csv", row.names = FALSE)

# Filter rows where adj.P.Val is less than or equal to 0.05
filtered_data <- combined_data %>%
  filter(Adjusted.P.value <= 0.05)
# filtered total keggs
write.csv(filtered_data, file = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_keggs_filtered.csv", row.names = FALSE)
filtered_data$Cell_type <- filtered_data$Directory
write.csv(table(filtered_data$Cell_type, filtered_data$DEG_test), file = "/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/soupx/PEBBLES_num_sig_KEGGs.csv", row.names = TRUE)
###########################################
## Create stacked bargraphs for the DEGs ##
###########################################
#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

#Create stacked bargraphs for the DEGs
counts <- as.data.frame(table(filtered_data$Directory, filtered_data$DEG_test))
# Assuming 'long_df' contains the relevant data
cluster_order <- c("L2_3_IT", "L4", "L5", "L6", "Pvalb", "Vip", "Sst", "Sncg", "Lamp5", "Peri", "Endo", "Oligo", "Astro", "Non-neuronal")

# Reorder the 'Cluster' factor levels
counts$Var1 <- factor(counts$Var1, levels = cluster_order)

ggplot(counts, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Number of KEGG pathways (adj.P.Val < 0.05)", title = "Number of KEGG Terms") +
  theme_minimal()+
  scale_fill_manual(values = polychrome_palette)+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    
    axis.text.x = element_text(angle = 40, size = 14, face = 'plain', hjust = 1.0, vjust = 1, family="sans", colour = 'black'),
    axis.text.y = element_text(angle = 0, size = 12, face = 'plain', vjust = 0.5, family="sans", colour = 'black'),
    axis.title = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.title.x = element_text(size = 14, face = 'plain', family="sans", colour = 'black', vjust = 1),
    axis.title.y = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "plain", family="sans"), # Text size
    title = element_text(size = 14, face = "plain", family="sans"))
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/05_Plots/num_sig_KEGG_terms.pdf",
                device = NULL,
                height = 10,
                width = 12)

# Filter KEGG terms
new_dataframe <- filtered_data %>%
  group_by(Term) %>%
  filter(
    any(DEG_test == "VehicleVsPCB_in_HETs") &
      any(DEG_test == "VehicleVsPCB_in_WTs") &
      any(DEG_test == "WTvsHET_in_PCBs") &
      any(DEG_test == "WTvsHET_in_Vehicle")
  ) %>%
  distinct()

# show the excluded terms
excluded_dataframe <- filtered_data %>%
  distinct() %>%
  anti_join(new_dataframe, by = "Term")

# Assuming filtered_data is your dataframe

# Get unique categories in DEG_test
categories <- unique(filtered_data$DEG_test)

# Create an empty list to store data frames
list_of_dataframes <- list()

# Iterate over categories
for (category in categories) {
  # Filter data for each category
  filtered_df <- filtered_data[filtered_data$DEG_test == category, ]
  # Append filtered data to the list
  list_of_dataframes[[category]] <- filtered_df
}

#Create input lists and make upset PDF
listInput <- list(VehicleVsPCB_in_HETs = list_of_dataframes$VehicleVsPCB_in_HETs$Term, 
                  VehicleVsPCB_in_WTs = list_of_dataframes$VehicleVsPCB_in_WTs$Term, 
                  WTvsHET_in_PCBs = list_of_dataframes$WTvsHET_in_PCBs$Term,
                  WTvsHET_in_Vehicle = list_of_dataframes$WTvsHET_in_Vehicle$Term)
pdf(glue("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/05_Plots/upset_pebbles_kegg.pdf"))
upset(fromList(listInput), sets = c('VehicleVsPCB_in_HETs','VehicleVsPCB_in_WTs', 'WTvsHET_in_PCBs','WTvsHET_in_Vehicle'), keep.order = TRUE)
dev.off()

ggVennDiagram(listInput) + scale_fill_gradient(low = 'grey', high = "red")
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/05_Plots/venn_pebbles_kegg.pdf",
                device = NULL,
                height = 10,
                width = 12)
ggVennDiagram(listInput, force_upset = TRUE, label = "count")

# Extract the "Term" column from the first dataframe in the list
common_terms <- as.character(list_of_dataframes[[1]]$Term)

# Iterate over the remaining dataframes in the list and find the intersection of "Term" values
for (i in 2:length(list_of_dataframes)) {
  common_terms <- intersect(common_terms, as.character(list_of_dataframes[[i]]$Term))
}

# Create an empty dataframe to store the intersecting rows
intersecting_df <- data.frame()

# Iterate over the list of dataframes and filter rows where "Term" matches the common terms
for (df in list_of_dataframes) {
  intersecting_df <- rbind(intersecting_df, df[df$Term %in% common_terms, ])
}

# Optionally, you may want to remove duplicate rows
intersecting_df <- unique(intersecting_df)
ggplot(intersecting_df, aes(x = Cell_type, y = Term, size = Odds.Ratio, color = Adjusted.P.value)) +
  geom_point() +
  scale_size_continuous(name = "Odds Ratio") +
  scale_color_continuous(name = "Adjusted P-Value") +
  labs(x = "Cell Type", y = "Term", title = "Dot Plot of Intersection Data") +
  theme_minimal() +
  scale_fill_manual(values = polychrome_palette)+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    axis.text.x = element_text(angle = 40, size = 14, face = 'plain', hjust = 1.0, vjust = 1, family="sans", colour = 'black'),
    axis.text.y = element_text(angle = 0, size = 12, face = 'plain', vjust = 0.5, family="sans", colour = 'black'),
    axis.title = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.title.x = element_text(size = 14, face = 'plain', family="sans", colour = 'black', vjust = 1),
    axis.title.y = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "plain", family="sans"), # Text size
    title = element_text(size = 14, face = "plain", family="sans"))
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/05_Plots/sig_KEGG_in_common.pdf",
                device = NULL,
                height = 10,
                width = 12)


# Strings common to all lists
common_keggs <- Reduce(intersect, listInput)

# Print common strings
print("Common kegg:")
print(common_keggs)

# Initialize a list to store unique terms for each list
unique_terms_list <- list()

# Find unique terms for each list and store them in the list
for (i in seq_along(listInput)) {
  unique_terms_list[[i]] <- setdiff(listInput[[i]], common_keggs)
}

# Print unique terms for each list
for (i in seq_along(unique_terms_list)) {
  cat("Unique terms in list", names(listInput)[i], ":", unique_terms_list[[i]], "\n")
}


matched_rows <- combined_data %>%
  filter(Term %in% unlist(unique_terms_list))
ggplot(matched_rows, aes(x = Directory, y = Term, size = Odds.Ratio, color = Adjusted.P.value)) +
  geom_point() +
  scale_size_continuous(name = "Odds Ratio") +
  scale_color_continuous(name = "Adjusted P-Value") +
  labs(x = "Cell Type", y = "Term", title = "Dot Plot of Intersection Data") +
  theme_minimal() +
  scale_fill_manual(values = polychrome_palette)+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    axis.text.x = element_text(angle = 40, size = 14, face = 'plain', hjust = 1.0, vjust = 1, family="sans", colour = 'black'),
    axis.text.y = element_text(angle = 0, size = 12, face = 'plain', vjust = 0.5, family="sans", colour = 'black'),
    axis.title = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.title.x = element_text(size = 14, face = 'plain', family="sans", colour = 'black', vjust = 1),
    axis.title.y = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "plain", family="sans"), # Text size
    title = element_text(size = 14, face = "plain", family="sans"))
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/05_Plots/sig_KEGG_in_common.pdf",
                device = NULL,
                height = 10,
                width = 12)

# Assuming 'filtered_data' is your dataframe
unique_strings <- unique(filtered_data$DEG_test)

# List to store individual dataframes
new_dataframes <- list()

# Loop through each unique string
for (string in unique_strings) {
  # Create a new dataframe for the current string
  new_df <- filtered_data[filtered_data$DEG_test == string, ]
  
  # Store the new dataframe in the list
  new_dataframes[[string]] <- new_df
}

# Now, you can access each dataframe using its corresponding unique string

top_5_rows <- filtered_data %>%
  filter(!is.na(Adjusted.P.value)) %>%  # Remove rows with NA in Adjusted.P.value
  arrange(desc(Odds.Ratio)) %>%         # Arrange rows by Adjusted.P.value
  group_by(DEG_test, Cell_type) %>%     # Group by DEG_test and Cell_type
  slice_head(n = 5)                     # Select the top 5 rows within each group

# View the top 5 rows based on Adjusted.P.value for each DEG_test and Cell_type
print(top_5_rows)

ggplot(table_result_sorted, aes(x = Cell_type, y = Term, size = Odds.Ratio, color = Adjusted.P.value)) +
  geom_point() +
  scale_size_continuous(name = "Odds Ratio") +
  scale_color_continuous(name = "Adjusted P-Value") +
  labs(x = "Cell Type", y = "Term", title = "Top 5 KEGG pathways for each condition") +
  theme_minimal() +
  scale_fill_manual(values = polychrome_palette)+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="sans", colour = 'black'),
    axis.text.x = element_text(angle = 40, size = 14, face = 'plain', hjust = 1.0, vjust = 1, family="sans", colour = 'black'),
    axis.text.y = element_text(angle = 0, size = 12, face = 'plain', vjust = 0.5, family="sans", colour = 'black'),
    axis.title = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.title.x = element_text(size = 14, face = 'plain', family="sans", colour = 'black', vjust = 1),
    axis.title.y = element_text(size = 14, face = 'plain', family="sans", colour = 'black'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "plain", family="sans"), # Text size
    title = element_text(size = 14, face = "plain", family="sans"))
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/05_Plots/sig_KEGG_in_common.pdf",
                device = NULL,
                height = 10,
                width = 12)

table_result <- table(top_5_rows$Term, top_5_rows$DEG_test)
# Convert the DEG_test column to factor preserving the order of levels
top_5_rows$DEG_test <- factor(top_5_rows$DEG_test, levels = unique(top_5_rows$DEG_test))

# Create the table with margins
table_result <- addmargins(table(top_5_rows$Term, top_5_rows$DEG_test))

# View the table
print(table_result)
ggplot2::ggsave("/Users/osman/Documents/GitHub/PEBBLES_mouse_snRNAseq/05_Plots/top5_sig_KEGG_perTest.pdf",
                device = NULL,
                height = 10,
                width = 12)
