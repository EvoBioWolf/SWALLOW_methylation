# Load necessary libraries
# If you don't have these installed, run:
# install.packages("ggplot2")
# install.packages("dplyr")
library(ggplot2)
library(dplyr)

# --- 1. Load the data ---
# IMPORTANT: This assumes your file is tab-separated. If it's comma-separated, change sep = ",".
# For very large files, consider using data.table::fread for faster loading:
# install.packages("data.table")
# library(data.table)
# df <- fread("your_large_data_file.txt")
# Otherwise, use read.delim:
df_all <- read.delim(gzfile("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/HR.85.maf2.ldwin100k.3500kb.inversion.ld.filter.gz"), header = T, sep = "", check.names = FALSE)
colnames(df_all) <- c("CHR_A",	"BP_A",	"SNP_A",	"CHR_B",	"BP_B",	"SNP_B",	"R2")

# --- 2. Data Preprocessing ---

# Select only the essential columns: SNP_A, BP_A, SNP_B, BP_B, and R2.
# We also filter the data to create the triangular heatmap.
# By keeping only rows where BP_A is less than or equal to BP_B,
# we ensure that each pair (SNP_A, SNP_B) is plotted only once,
# forming the lower triangle of the heatmap.
df_filtered <- df_all %>%
  dplyr::select(SNP_A, BP_A, SNP_B, BP_B, R2) %>%
  dplyr::filter(BP_A <= BP_B) # This creates the lower triangular heatmap

# To ensure correct ordering of SNPs on the axes (which is crucial for genomic heatmaps),
# we need to convert SNP_A and SNP_B into factors and order their levels
# based on their base pair positions (BP_A and BP_B).
# First, create a unique list of all SNPs and their positions, then sort by position.
all_snps_positions <- df_filtered %>%
  dplyr::select(SNP_A, BP_A) %>%
  dplyr::rename(SNP = SNP_A, BP = BP_A) %>% # Rename for consistent column names
  bind_rows( # Combine with SNPs from the B column
    df_filtered %>%
      dplyr::select(SNP_B, BP_B) %>%
      dplyr::rename(SNP = SNP_B, BP = BP_B)
  ) %>%
  distinct(SNP, BP) %>% # Get unique SNP-BP pairs
  arrange(BP) # Order by base pair position

# Now, apply this ordered list of SNPs as factor levels to both SNP_A and SNP_B.
df_filtered$SNP_A <- factor(df_filtered$SNP_A, levels = all_snps_positions$SNP)
df_filtered$SNP_B <- factor(df_filtered$SNP_B, levels = all_snps_positions$SNP)


df_filtered_high <- df_filtered %>%
  filter(R2 > 0.6)

df_filtered_low <- df_filtered %>%
  filter(R2 <0.6)

df_filtered_low <- df_filtered_low %>%
  sample_frac(0.20)

# --- 3. Plotting with ggplot2 ---

# Define a custom color palette that goes from blue to orange, similar to your example.
# You can adjust these hex codes or use RColorBrewer palettes for more options.
ld_colors <- c("#ADD8E6", "#4682B4", "#104E8B", "#FFA500", "#FF4500", "#8B0000") # Light blue to dark red/orange

p <- ggplot(df_filtered, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  # Use geom_tile to create the heatmap cells.
  # 'color = "white"' adds thin white borders between tiles for better visual separation.
  geom_tile() +
  # Define the fill color scale for R2 values.
  # scale_fill_gradientn allows for multiple colors in the gradient.
  # limits = c(0, 1) sets the range for R2 values.
  # breaks defines where the labels on the color key will appear.
  # name sets the title for the color legend, using expression for R^2 formatting.
  scale_fill_gradientn(colors = ld_colors,
                       limits = c(0, 1),
                       breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
                       name = expression(R^2 ~ "Color Key")) +
  # Use a minimal theme for a clean background.
  theme_minimal() +
  # Customize labels and titles.
  labs(title = "LD Heatmap",
       x = NULL, # Remove x-axis label as SNP names will serve as labels
       y = NULL) + # Remove y-axis label
  theme(
    # Customize axis text appearance and rotation for readability.
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    # Center the plot title and make it bold.
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    # Customize the legend position and appearance.
    legend.position = "bottom", # Place the legend below the plot
    legend.key.width = unit(2, "cm"), # Make the legend bar wider
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    # Remove grid lines for a cleaner heatmap look.
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove the background panel and plot background.
    panel.background = element_blank(),
    plot.background = element_blank(),
    # Adjust plot margins if needed.
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Top, Right, Bottom, Left margins
  ) +
  # Ensure the tiles are square by fixing the aspect ratio.
  coord_fixed(ratio = 1)

# Print the plot to display it.
print(p)

# Optional: Save the plot to a file (e.g., PNG).
# Adjust width, height, and dpi as needed for desired output quality.
ggsave("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/2_all_subset_snp_chr11_inv_heatmap.png", plot = p, width = 8, height = 8, dpi = 600)



###PROCESS TO PLOT SNPs and Meth separately
# --- 2. Data Preprocessing ---

#For proper plotting, where I will order the heatmap by SNPS and DNAm loci, I need to make sure all methylation - snp comparisons are in the same order
#The following filtering ensures all S_ values are in A columns and all M_ values are in B columns
# Perform the conditional swap using ifelse for each column
df_all_swapped <- df_all %>%
  dplyr::mutate(
    # Define the condition for swapping: SNP_A starts with "M_" AND SNP_B starts with "S_"
    condition_to_swap = startsWith(SNP_A, "M_") & startsWith(SNP_B, "S_"),
    
    # Use case_when to conditionally assign values for SNP_A
    # If condition_to_swap is TRUE, the new SNP_A takes the original SNP_B value.
    # Otherwise, it keeps its original SNP_A value.
    new_SNP_A = case_when(
      condition_to_swap ~ SNP_B,
      TRUE ~ SNP_A # Default case: keep original SNP_A
    ),
    # Conditionally assign values for BP_A, mirroring the SNP_A logic
    new_BP_A = case_when(
      condition_to_swap ~ BP_B,
      TRUE ~ BP_A
    ),
    # Conditionally assign values for SNP_B
    # If condition_to_swap is TRUE, the new SNP_B takes the original SNP_A value.
    # Otherwise, it keeps its original SNP_B value.
    new_SNP_B = case_when(
      condition_to_swap ~ SNP_A, # This SNP_A refers to the original SNP_A from df_all
      TRUE ~ SNP_B
    ),
    # Conditionally assign values for BP_B, mirroring the SNP_B logic
    new_BP_B = case_when(
      condition_to_swap ~ BP_A, # This BP_A refers to the original BP_A from df_all
      TRUE ~ BP_B
    )
  ) %>%
  # Select the newly created columns as the main SNP_A, BP_A, SNP_B, BP_B
  # and retain the R2 column, discarding the temporary 'condition_to_swap' column.
  dplyr::select(SNP_A = new_SNP_A, BP_A = new_BP_A,
                SNP_B = new_SNP_B, BP_B = new_BP_B,
                R2)

# Select only the essential columns: SNP_A, BP_A, SNP_B, BP_B, and R2.
# We also filter the data to create the triangular heatmap.
# By keeping only rows where BP_A is less than or equal to BP_B,
# we ensure that each pair (SNP_A, SNP_B) is plotted only once,
# forming the lower triangle of the heatmap.
df_filtered <- df_all_swapped %>%
  dplyr::select(SNP_A, BP_A, SNP_B, BP_B, R2) 

# To ensure correct ordering of SNPs on the axes (which is crucial for genomic heatmaps),
# we need to convert SNP_A and SNP_B into factors and order their levels
# based on their base pair positions (BP_A and BP_B).
# First, create a unique list of all SNPs and their positions from both SNP_A and SNP_B.
all_snps_positions <- df_filtered %>%
  dplyr::select(SNP_A, BP_A) %>%
  dplyr::rename(SNP = SNP_A, BP = BP_A) %>% # Rename for consistent column names
  bind_rows( # Combine with SNPs from the B column
    df_filtered %>%
      dplyr::select(SNP_B, BP_B) %>%
      dplyr::rename(SNP = SNP_B, BP = BP_B)
  ) %>%
  distinct(SNP, BP) # Get unique SNP-BP pairs

# Separate 'S_' and 'M_' SNPs and order them by base pair position
s_snps_positions <- all_snps_positions %>%
  filter(startsWith(SNP, "S_")) %>%
  arrange(BP)

m_snps_positions <- all_snps_positions %>%
  filter(startsWith(SNP, "M_")) %>%
  arrange(BP)

# Combine the ordered 'S_' SNPs and 'M_' SNPs.
# The 'S_' SNPs will come first, followed by 'M_' SNPs.
ordered_all_snps <- c(s_snps_positions$SNP, m_snps_positions$SNP)

# Now, apply this custom ordered list of SNPs as factor levels to both SNP_A and SNP_B.
df_filtered$SNP_A <- factor(df_filtered$SNP_A, levels = ordered_all_snps)
df_filtered$SNP_B <- factor(df_filtered$SNP_B, levels = ordered_all_snps)

# Determine the position for the visual break between S_ and M_ markers.
# This will be just after the last S_ SNP in the ordered list.
s_snp_count <- nrow(s_snps_positions)

# --- 3. Plotting with ggplot2 ---

# Define a custom color palette that goes from blue to orange, similar to your example.
ld_colors <- c("#ADD8E6", "#4682B4", "#104E8B", "#FFA500", "#FF4500", "#8B0000") # Light blue to dark red/orange
ld_colors <- c("#FFFFE0", "#FFDDAA", "#FFB03B", "#FFA500", "#FF4500", "#8B0000")

p <- ggplot(df_filtered, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  # Use geom_tile to create the heatmap cells.
  geom_tile() +
  # Add vertical and horizontal lines to create a visual break between S_ and M_ markers.
  # The lines are placed at s_snp_count + 0.5 to appear between the last S_ tile and the first M_ tile.
  geom_vline(xintercept = s_snp_count + 0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = s_snp_count + 0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  # Define the fill color scale for R2 values.
  scale_fill_gradientn(colors = ld_colors,
                       limits = c(0, 1),
                       breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
                       name = expression(R^2 ~ "Color Key")) +
  # Use a minimal theme for a clean background.
  theme_minimal() +
  # Customize labels and titles.
  labs(title = "LD Heatmap (S_ vs M_ Sorted)", # Updated title for clarity
       x = "SNPs", # Generic x-axis label
       y = "SNPs") + # Generic y-axis label
  theme(
    # Customize axis text appearance and rotation for readability.
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    # Center the plot title and make it bold.
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    # Customize the legend position and appearance.
    legend.position = "bottom", # Place the legend below the plot
    legend.key.width = unit(2, "cm"), # Make the legend bar wider
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    # Remove grid lines for a cleaner heatmap look.
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove the background panel and plot background.
    panel.background = element_blank(),
    plot.background = element_blank(),
    # Adjust plot margins if needed.
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Top, Right, Bottom, Left margins
  ) +
  # Ensure the tiles are square by fixing the aspect ratio.
  coord_fixed(ratio = 1)

# Print the plot to display it.
print(p)

# Optional: Save the plot to a file (e.g., PNG).
# Adjust width, height, and dpi as needed for desired output quality.
ggsave("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/2_all_ordered_chr11_inv_heatmap.png", plot = p, width = 20, height = 20, dpi = 600)

# --- New code to filter based on R2 values for SNP_A and SNP_B pairs ---

# Identify SNP_A and SNP_B pairs where at least one R2 value is greater than 0.6
snps_to_keep_based_on_R2 <- df_filtered %>%
  group_by(SNP_A, SNP_B) %>%
  filter(any(R2 > 0.6)) %>%
  ungroup() %>%
  select(SNP_A, SNP_B) %>%
  distinct()

# Filter df_filtered_modified to keep only rows corresponding to these identified pairs
df_filtered_modified <- df_filtered %>%
  inner_join(snps_to_keep_based_on_R2, by = c("SNP_A", "SNP_B"))

df_filtered_modified$SNP_A <- as.character(df_filtered_modified$SNP_A)
df_filtered_modified$SNP_B <- as.character(df_filtered_modified$SNP_B)

##get M_ and S_ in the right positions
df_all_swapped <- df_filtered_modified %>%
  dplyr::mutate(
    # Define the condition for swapping: SNP_A starts with "M_" AND SNP_B starts with "S_"
    condition_to_swap = startsWith(SNP_A, "M_") & startsWith(SNP_B, "S_"),
    
    # Use case_when to conditionally assign values for SNP_A
    # If condition_to_swap is TRUE, the new SNP_A takes the original SNP_B value.
    # Otherwise, it keeps its original SNP_A value.
    new_SNP_A = case_when(
      condition_to_swap ~ SNP_B,
      TRUE ~ SNP_A # Default case: keep original SNP_A
    ),
    # Conditionally assign values for BP_A, mirroring the SNP_A logic
    new_BP_A = case_when(
      condition_to_swap ~ BP_B,
      TRUE ~ BP_A
    ),
    # Conditionally assign values for SNP_B
    # If condition_to_swap is TRUE, the new SNP_B takes the original SNP_A value.
    # Otherwise, it keeps its original SNP_B value.
    new_SNP_B = case_when(
      condition_to_swap ~ SNP_A, # This SNP_A refers to the original SNP_A from df_all
      TRUE ~ SNP_B
    ),
    # Conditionally assign values for BP_B, mirroring the SNP_B logic
    new_BP_B = case_when(
      condition_to_swap ~ BP_A, # This BP_A refers to the original BP_A from df_all
      TRUE ~ BP_B
    )
  ) %>%
  # Select the newly created columns as the main SNP_A, BP_A, SNP_B, BP_B
  # and retain the R2 column, discarding the temporary 'condition_to_swap' column.
  dplyr::select(SNP_A = new_SNP_A, BP_A = new_BP_A,
                SNP_B = new_SNP_B, BP_B = new_BP_B,
                R2)


# Select only the essential columns: SNP_A, BP_A, SNP_B, BP_B, and R2.
# We also filter the data to create the triangular heatmap.
# By keeping only rows where BP_A is less than or equal to BP_B,
# we ensure that each pair (SNP_A, SNP_B) is plotted only once,
# forming the lower triangle of the heatmap.
df_filtered_modified <- df_all_swapped %>%
  dplyr::select(SNP_A, BP_A, SNP_B, BP_B, R2) 

# To ensure correct ordering of SNPs on the axes (which is crucial for genomic heatmaps),
# we need to convert SNP_A and SNP_B into factors and order their levels
# based on their base pair positions (BP_A and BP_B).
# First, create a unique list of all SNPs and their positions from both SNP_A and SNP_B.
all_snps_positions <- df_filtered_modified %>%
  dplyr::select(SNP_A, BP_A) %>%
  dplyr::rename(SNP = SNP_A, BP = BP_A) %>% # Rename for consistent column names
  bind_rows( # Combine with SNPs from the B column
    df_filtered_modified %>%
      dplyr::select(SNP_B, BP_B) %>%
      dplyr::rename(SNP = SNP_B, BP = BP_B)
  ) %>%
  distinct(SNP, BP) # Get unique SNP-BP pairs

# Separate 'S_' and 'M_' SNPs and order them by base pair position
s_snps_positions <- all_snps_positions %>%
  filter(startsWith(SNP, "S_")) %>%
  arrange(BP)

m_snps_positions <- all_snps_positions %>%
  filter(startsWith(SNP, "M_")) %>%
  arrange(BP)

# Combine the ordered 'S_' SNPs and 'M_' SNPs.
# The 'S_' SNPs will come first, followed by 'M_' SNPs.
ordered_all_snps <- c(s_snps_positions$SNP, m_snps_positions$SNP)

# Now, apply this custom ordered list of SNPs as factor levels to both SNP_A and SNP_B.
df_filtered$SNP_A <- factor(df_filtered$SNP_A, levels = ordered_all_snps)
df_filtered$SNP_B <- factor(df_filtered$SNP_B, levels = ordered_all_snps)

# Determine the position for the visual break between S_ and M_ markers.
# This will be just after the last S_ SNP in the ordered list.
s_snp_count <- nrow(s_snps_positions)

# --- 3. Plotting with ggplot2 ---

# Define a custom color palette that goes from blue to orange, similar to your example.
ld_colors <- c("#ADD8E6", "#4682B4", "#104E8B", "#FFA500", "#FF4500", "#8B0000") # Light blue to dark red/orange
ld_colors <- c("#FFFFE0", "#FFDDAA", "#FFB03B", "#FFA500", "#FF4500", "#8B0000")

p <- ggplot(df_filtered_modified, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  # Use geom_tile to create the heatmap cells.
  geom_tile() +
  # Add vertical and horizontal lines to create a visual break between S_ and M_ markers.
  # The lines are placed at s_snp_count + 0.5 to appear between the last S_ tile and the first M_ tile.
  geom_vline(xintercept = s_snp_count + 0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = s_snp_count + 0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  # Define the fill color scale for R2 values.
  scale_fill_gradientn(colors = ld_colors,
                       limits = c(0, 1),
                       breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
                       name = expression(R^2 ~ "Color Key")) +
  # Use a minimal theme for a clean background.
  theme_minimal() +
  # Customize labels and titles.
  labs(title = "LD Heatmap (S_ vs M_ Sorted)", # Updated title for clarity
       x = "SNPs", # Generic x-axis label
       y = "SNPs") + # Generic y-axis label
  theme(
    # Customize axis text appearance and rotation for readability.
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    # Center the plot title and make it bold.
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    # Customize the legend position and appearance.
    legend.position = "bottom", # Place the legend below the plot
    legend.key.width = unit(2, "cm"), # Make the legend bar wider
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    # Remove grid lines for a cleaner heatmap look.
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove the background panel and plot background.
    panel.background = element_blank(),
    plot.background = element_blank(),
    # Adjust plot margins if needed.
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Top, Right, Bottom, Left margins
  ) +
  # Ensure the tiles are square by fixing the aspect ratio.
  coord_fixed(ratio = 1)

# Print the plot to display it.
print(p)

# Optional: Save the plot to a file (e.g., PNG).
# Adjust width, height, and dpi as needed for desired output quality.
ggsave("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/2_all_ordered_chr11_inv_heatmap.png", plot = p, width = 20, height = 20, dpi = 600)
