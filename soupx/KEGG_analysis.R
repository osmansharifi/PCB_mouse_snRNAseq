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
