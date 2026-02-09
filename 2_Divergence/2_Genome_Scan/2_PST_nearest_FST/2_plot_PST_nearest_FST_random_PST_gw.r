.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))

library(vegan)
library(ggplot2, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
library(plyr)
library(broom)
library(readxl)
library(cluster)
library(ape)
library(languageR)
library(packfor)
library(ggrepel)
library(knitr)
library(ggvenn)
library(ggbreak)
library(gridExtra)
library(tidyr)
library(dplyr)
library(purrr)

###Plotting from observed and simulated mean distances in CpGislands, shores, shelves
## Currently, plotting from ~36,000 regions

####SET-UP###########

#setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_site_pst/")

#Load in PST data - it was either run on the top 1% or top 0.1% of PST outliers
pst_data <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_site_pst/0_pst_variant_c05.txt", header=T)

# Reshape the PST data into long form
pst_data_long <- pst_data %>%
  pivot_longer(-site.id, names_to = "comparison", values_to = "pst") %>%
  drop_na() %>%
  separate(comparison, into = c("pop1", "pop2"), sep = "\\.vs\\.") %>%
  mutate(pop2 = gsub("\\.pst$", "", pop2))


pst_data_long <- pst_data_long %>%
  separate(site.id, into = c("chr", "pos", "end"), sep = "_") %>%
  mutate(pos = as.numeric(pos))

#pop key function
generate_pop_key <- function(pop1, pop2) {
  sorted_pops <- sort(c(as.character(pop1), as.character(pop2)))
  return(paste(sorted_pops[1], sorted_pops[2], sep = "_"))
}

# Identify top 1% PST loci for each combination of pop1 and pop2
top_1_pct_pst <- pst_data_long %>%
  group_by(pop1,pop2) %>%
  arrange(desc(pst)) %>%
  mutate(rank = row_number()) %>%
  group_by(pop1,pop2) %>%
  filter(pst >= quantile(pst, 0.99, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-rank) # Remove rank column if not needed

#add pop key
top_1_pct_pst <- top_1_pct_pst %>%
  rowwise() %>%
  mutate(pop_key = generate_pop_key(pop1, pop2)) %>%
  ungroup()

#write.table(top_1_pct_pst, file= "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/1_top_pst_99_results.txt", quote=F, row.names=F)
#top_1_pct_pst <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_simulation_results/1_top_pst_999_results.txt", header=T)

########Plotting where CpGi show up in genome, how many outliers per chromosome, etc. 
# Count the occurrences of each chromosome
chr_counts <- top_1_pct_pst %>%
  group_by(chr) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  # Optionally, order chromosomes for better visualization (e.g., numerically for chr)
  mutate(chr = factor(chr, levels = mixedsort(unique(chr))))

# Create the bar plot
chromosome_barplot <- ggplot(chr_counts, aes(x = chr, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") + # Use "identity" for pre-counted data
  labs(
    title = "Frequency of Top 0.1% PST Loci by Chromosome",
    x = "Chromosome",
    y = "Number of Loci"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5) # Center the title
  )

# Print the plot
print(chromosome_barplot)

##Add genome information
hr.genome <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome")
colnames(hr.genome) <- c("chr", "max_bp")
hr.genome <- hr.genome %>%
  filter(!grepl("_random|W|Z", chr))


##Correct for chromosome size
chr_counts_with_length <- left_join(chr_counts, hr.genome, by = "chr") %>%
  # Ensure all necessary chromosomes are present
  filter(!is.na(max_bp))

# 4. Calculate the normalized count (loci per megabase)
normalized_chr_counts <- chr_counts_with_length %>%
  mutate(normalized_count = count / (max_bp / 1e6)) %>%
  # Order chromosomes for better plotting (e.g., numerically)
  mutate(chr = factor(chr, levels = gtools::mixedsort(unique(chr))))

# 5. Create the bar plot with the normalized values
normalized_barplot <- ggplot(normalized_chr_counts, aes(x = chr, y = normalized_count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(
    title = "Normalized Frequency of Top 0.1% PST Loci",
    subtitle = "Loci per Megabase of Chromosome Length",
    x = "Chromosome",
    y = "Normalized Count (Loci per Mb)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

# Print the final normalized plot
print(normalized_barplot)

#########SET-UP SIMULATION RESULTS###########

final_results_gw <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_PST_nearest_FST/2_simulation_results/1_PST_99_FST_99_1000sim_all_results.txt", header=T)

##Correct p-values based on multiple pairwise comparisons
final_results_gw$corrected_p_values <- p.adjust(p = final_results_gw$p_value, method = "bonferroni")

final_results_long <- final_results_gw %>%
  pivot_longer(
    cols = starts_with("Sim_"),  # Pivot all simulation columns
    names_to = "simulation",    # New column for simulation identifiers
    values_to = "simulated_value"  # New column for simulation values
  )

# Plot histogram for each combination of pop_key and chr
ggplot(final_results_long, aes(x = simulated_value)) +
  geom_histogram(bins = 30, fill = "skyblue", alpha = 0.6, color = "black") +
  geom_vline(
    data = final_results_long,
    aes(xintercept = observed_mean_distance, color = "Observed Mean"),
    linetype = "dashed",
    size = 1
  ) +
  facet_wrap(~ pop_key, scales = "free") +  # Facet by pop_key and chr
  theme_minimal() +
  labs(
    title = "Distribution of Simulation Values with Observed Mean",
    x = "Normalized Distance",
    y = "Frequency"
  ) +
  scale_color_manual(values = "grey33") +
  theme(legend.position = "bottom")

###Check number of pairwise comparisons that have distances below the 5% or 1% threshold
threshold_results <- final_results_long %>%
  group_by(pop_key) %>%
  summarise(
    threshold_5 = quantile(simulated_value, probs = 0.05, na.rm = TRUE),
    threshold_1 = quantile(simulated_value, probs = 0.01, na.rm = TRUE),
    threshold_95 = quantile(simulated_value, probs = 0.95, na.rm = TRUE),
    observed_mean_distance = unique(observed_mean_distance)
  ) %>%
  mutate(
    below_5_percent = observed_mean_distance < threshold_5,
    below_1_percent = observed_mean_distance < threshold_1,
    above_95_percent = observed_mean_distance > threshold_95
  )

# Calculate proportion of TRUE for below_5_percent
prop_below_5 <- sum(threshold_results$below_5_percent, na.rm = TRUE) / nrow(threshold_results)
prop_below_1 <- sum(threshold_results$below_1_percent, na.rm = TRUE) / nrow(threshold_results)

# Calculate proportion of TRUE for above_95_percent
prop_above_95 <- sum(threshold_results$above_95_percent, na.rm = TRUE) / nrow(threshold_results)

# Display results
cat("Proportion of observed means below the 5% threshold:", prop_below_5, "\n")
#57% of comparisons are below 5% threshold
cat("Proportion of observed means below the 1% threshold:", prop_below_1, "\n")
#42% ARE BELOW THE 1% THRESHOLD
cat("Proportion of observed means above the 95% threshold:", prop_above_95, "\n")
#9% are above 95% threshold

##SUmmaries of p-value and adjusted p-values
summary(final_results_gw$corrected_p_values)
summary(final_results_gw$p_value)


###ADD in other comparisons to look

# 1. Add in data on strength of FST based on whole genome estimates 
fst.gw <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/1_Population_Structure/4_pixy/2_FST/1_fst_avg_5kb_gw.txt", header=T)
fst.gw <- fst.gw %>%
  rowwise() %>%
  dplyr::mutate(pop_key = generate_pop_key(pop1, pop2)) %>%
  ungroup()
final_results_gw <- left_join(final_results_gw, fst.gw, by=c("pop_key"))

# 2.Create matrix of observed mean and fst (originally for Mantel tests...not sure if needed now)
obs_mean <- final_results_gw$observed_mean_distance
fst <- final_results_gw$avg_fst

# 3. Calculate the correlation
correlation.fst <- cor(obs_mean, fst)
cor.test(obs_mean, fst)

# 4. Create the scatterplot for FST
ggplot(final_results_gw, aes(x = avg_fst, y = observed_mean_distance)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) + # Add regression line
  theme_minimal() +
  labs(
    title = paste("Correlation between Mean distance to FST peaks and genomewide FST"),
    subtitle = paste("Correlation:", round(correlation.fst, 3)),
    x = "Pairwise FST",
    y = "Observed Mean Distance"
  )

# 5. Order the pop_key based on FST
ordered_pop_keys <- final_results_gw %>%
  arrange(fst) %>%
  pull(pop_key)

# 6. Prepare the data
final_results_plot <- final_results_gw %>%
  mutate(pop_key = factor(pop_key, levels = ordered_pop_keys))

# 4. Add this in a way that can be plotted
final_results_plot <- final_results_plot %>%
  mutate(pop_key = factor(pop_key, levels = ordered_pop_keys), # Make sure pop_key is a factor and ordered
         fst_category = case_when(
           fst >= 0 & fst < 0.01 ~ "None",
           fst >= 0.01 & fst < 0.02 ~ "Low",
           fst >= 0.02 & fst < 0.04 ~ "Medium",
           fst >= 0.03 & fst <= 0.06 ~ "High",  # Changed to be inclusive of 0.4
           TRUE ~ NA_character_ #handles NA and values outside of the specified ranges.
         ),
         fst_category = factor(fst_category, levels = c("None", "Low", "Medium", "High")) #important to maintain order
  )

final_results_plot_sub <- final_results_plot %>%
  filter(pop_key %in% c("gutturalis_savignii", "gutturalis_rustica.r", "rustica.r_tytleri", "savignii_transitiva")) %>%
  mutate(plot_number = c(4, 3, 2, 1))

# Define the color palette.  This uses a named vector, which is good practice.
fst_colors <- c("None" = "grey88", "Low" = "grey66", "Medium" = "grey11", "High" = "grey11")


#Get labels ready for plotting
custom_x_breaks <- c("gutturalis_savignii", "gutturalis_rustica.r", "rustica.r_tytleri", "savignii_transitiva")
custom_x_labels <- c("GU - SA", "GU - RUr", "RUr - TY", "SA - TR") # Your desired display names

# 4. Create the boxplot with the localized horizontal line
l <- ggplot(final_results_plot_sub, aes(x = pop_key, y = observed_mean_distance)) +
  geom_boxplot(aes(y = Distance, colour= fst_category), fill = "grey99",
               data = final_results_plot_sub %>%
                 dplyr::select(pop_key, fst_category, starts_with("Sim_")) %>%
                 pivot_longer(cols = starts_with("Sim_"),
                              names_to = "Simulation",
                              values_to = "Distance")) +
  scale_colour_manual(values = fst_colors, na.value = "grey") +
  geom_segment(aes(x = as.numeric(plot_number) - 0.4,  # Adjust for line start position
                   xend = as.numeric(plot_number) + 0.4, # Adjust for line end position
                   y = observed_mean_distance,
                   yend = observed_mean_distance),
               color = "orange", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = 4.5,
                ymin = 0.195, ymax = 0.21),
            fill = "white", color = NA, # No border
            inherit.aes = FALSE,
            data = NULL) + 
  geom_text(data = data.frame(pop_key = 4.0, label = "p-value"),
            aes(y = .20, x=4.3,
                label = label),
            inherit.aes = FALSE,       # Don't inherit main plot aesthetics
            hjust = 0.5, vjust = 0,    # Center horizontally, align to top vertically
            size = 4, fontface = "bold") +
  geom_text(aes(y = .20, 
                x = pop_key,            
                label = ifelse(corrected_p_values < 0.001, "<0.001", sprintf("%.3f", corrected_p_values))),
            hjust = 0.5, vjust = 0.5, # Center alignment for values
            size = 4) + 
  theme_minimal() +
  labs(
    y = "Distance to nearest FST outliers"
  ) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0.14, 0.18, 0.22)) +
  scale_x_discrete(breaks = custom_x_breaks, labels = custom_x_labels) + 
  coord_flip()

# Print the plot
print(l)

#save plot
#ggsave("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_PST_nearest_FST/3_Figure2b_simulation_results_subset.png", plot = l, width = 5, height = 2.5, dpi = 600)


#######Plot p-values for figure
final_results_plot <- final_results_plot %>%
  mutate(plot.p = "All Comparisons")

##Plot only p-values
p_value_plot <- ggplot(final_results_plot, aes(x = plot.p, y = corrected_p_values)) +
  geom_boxplot(color = "grey44", size = 1) + # Blue points for p-values
  theme_minimal() +
  labs(
    x = "Population Key",
    y = "P-value"
  ) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position = "none") +
  coord_flip() # Flip coordinates for horizontal orientation

# Print the plot
print(p_value_plot)

#save plot
#ggsave("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_simulation_results/3_fig_simulation_pvalues.png", plot = p_value_plot, width = 7, height = 1, dpi = 600)


#####PLOT ALL POPS - for supplements #######
#shorten names
final_results_plot <- final_results_plot %>%
  mutate(
    pop_key1 = str_split(pop_key, "_", simplify = TRUE) %>% 
      as.data.frame() %>%
      mutate(
        V1 = case_when(
          str_detect(V1, "^rustica\\.eu") ~ "RUe",
          str_detect(V1, "^rustica\\.r")  ~ "RUr",
          TRUE ~ toupper(substr(V1, 1, 2))
        ),
        V2 = case_when(
          str_detect(V2, "^rustica\\.eu") ~ "RUe",
          str_detect(V2, "^rustica\\.r")  ~ "RUr",
          TRUE ~ toupper(substr(V2, 1, 2))
        ),
        pop_key1 = paste0(V1,"_", V2)
      ) %>%
      pull(pop_key1)
  )

final_results_plot <- final_results_plot %>%
    mutate(pop_key1 = fct_reorder(pop_key1, avg_fst, .fun = max, .desc = FALSE))

# Define the color palette.  This uses a named vector, which is good practice.
fst_colors <- c("None" = "grey88", "Low" = "grey66", "Medium" = "grey44", "High" = "grey11")
#Create the boxplot with the localized horizontal line
p <- ggplot(final_results_plot, aes(x = pop_key1, y = observed_mean_distance)) +
  geom_boxplot(aes(y = Distance, colour= fst_category), fill = "grey88",
               data = final_results_plot %>%
                 dplyr::select(pop_key1, fst_category, avg_fst, starts_with("Sim_")) %>%
                 pivot_longer(cols = starts_with("Sim_"),
                              names_to = "Simulation",
                              values_to = "Distance")) +
  ylim(0,0.33)+
  scale_colour_manual(values = fst_colors, na.value = "grey") +
  geom_segment(aes(x = as.numeric(as.factor(final_results_plot$pop_key1)) - 0.4,  # Adjust for line start position
                   xend = as.numeric(as.factor(final_results_plot$pop_key1)) + 0.4, # Adjust for line end position
                   y = observed_mean_distance,
                   yend = observed_mean_distance),
               color = "orange", linewidth = 1) +
  theme_minimal() +
  labs(
    y = "Mean Distance to nearest FST outliers"
  ) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none") +
  coord_flip()

p


##PLot with p-values and FST values
p <- ggplot(final_results_plot, aes(x = pop_key1, y = observed_mean_distance)) +
  geom_boxplot(aes(y = Distance, colour= fst_category), fill = "grey88",
               data = final_results_plot %>%
                 dplyr::select(pop_key1, fst_category, starts_with("Sim_")) %>%
                 pivot_longer(cols = starts_with("Sim_"),
                              names_to = "Simulation",
                              values_to = "Distance")) +
  scale_colour_manual(values = fst_colors, na.value = "grey") +
  geom_segment(aes(x = as.numeric(as.factor(final_results_plot$pop_key1)) - 0.4,  # Adjust for line start position
                   xend = as.numeric(as.factor(final_results_plot$pop_key1)) + 0.4, # Adjust for line end position
                   y = observed_mean_distance,
                   yend = observed_mean_distance),
               color = "orange", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = 22,
                ymin = 0.24, ymax = 0.33),
            fill = "white", color = NA, # No border
            inherit.aes = FALSE,
            data = NULL) + 
  geom_text(data = data.frame(pop_key1 = 4.0, label = "p-value"),
            aes(y = .27, x=21.5,
                label = label),
            inherit.aes = FALSE,       # Don't inherit main plot aesthetics
            hjust = 0.5, vjust = 0,    # Center horizontally, align to top vertically
            size = 4, fontface = "bold") +
  geom_text(data = data.frame(pop_key1 = 4.0, label = "FST"),
            aes(y = .30, x=21.5,
                label = label),
            inherit.aes = FALSE,       # Don't inherit main plot aesthetics
            hjust = 0.5, vjust = 0,    # Center horizontally, align to top vertically
            size = 4, fontface = "bold") +
  geom_text(aes(y = .30, 
                x = pop_key1,            
                label = sprintf("%.3f", avg_fst)),
            hjust = 0.5, vjust = 0.5, # Center alignment for values
            size = 4) + 
  geom_text(aes(y = .27, 
                x = pop_key1,            
                label = ifelse(corrected_p_values < 0.001, "<0.001", sprintf("%.3f", corrected_p_values))),
            hjust = 0.5, vjust = 0.5, # Center alignment for values
            size = 4) + 
  theme_minimal() +
  labs(
    y = "Distance to nearest FST outliers"
  ) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0.14, 0.18, 0.22)) +
  #scale_x_discrete(breaks = custom_x_breaks, labels = custom_x_labels) + 
  coord_flip()

# Print the plot
print(p)

#ggsave("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_PST_nearest_FST/3_SuppFig_simulation_all.png", plot = p, width = 12, height = 6, dpi = 600)
