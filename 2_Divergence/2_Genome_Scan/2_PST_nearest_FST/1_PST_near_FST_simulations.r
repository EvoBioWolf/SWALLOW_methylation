
# --- Load libraries ---
.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))

library(vegan)
#library(ggplot2, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
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

# --- Working Directory ---
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_site_pst/")

# --- Input Data ---
pst_data <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_site_pst/0_pst_variant_c05.txt", header=T)
fst_data <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/1_Population_Structure/4_pixy/HR.168.invariant.10kb/HR.168.invariant.10kb_fst.txt", header = TRUE)
colnames(fst_data) <- c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2", "avg_wc_fst", "no_snps")

#genome information
hr.genome <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome")
colnames(hr.genome) <- c("chr", "max_bp")
hr.genome <- hr.genome %>%
  filter(!grepl("_random|W|Z", chr))

# --- Functions ---
generate_pop_key <- function(pop1, pop2) {
  sorted_pops <- sort(c(as.character(pop1), as.character(pop2)))
  return(paste(sorted_pops[1], sorted_pops[2], sep = "_"))
}

# --- Format PST data ---
#convert to long form so pst for each comparison is separate
pst_data_long <- pst_data %>%
  pivot_longer(-site.id, names_to = "comparison", values_to = "pst") %>%
  drop_na() %>%
  separate(comparison, into = c("pop1", "pop2"), sep = "\\.vs\\.") %>%
  mutate(pop2 = gsub("\\.pst$", "", pop2))

#separate site.id column 
pst_data_long <- pst_data_long %>%
  separate(site.id, into = c("chr", "pos", "end"), sep = "_") %>%
  mutate(pos = as.numeric(pos))

# Identify top 1% PST loci
top_1_pct_pst <- pst_data_long %>%
  group_by(pop1,pop2) %>%
  arrange(desc(pst)) %>%
  mutate(rank = row_number()) %>%
  group_by(pop1,pop2) %>%
  filter(pst >= quantile(pst, 0.99, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-rank) # Remove rank column if not needed

#add the Pop_key for comparison to FST
top_1_pct_pst <- top_1_pct_pst %>%
  rowwise() %>%
  mutate(pop_key = generate_pop_key(pop1, pop2)) %>%
  ungroup()

# --- Format PST data ---
#filter top 1% FST windows in each pop comparison
fst_peaks <- fst_data %>%
  group_by(pop1, pop2, chromosome) %>% # Group by pop1, pop2, and chromosome
  mutate(quantile_threshold = quantile(avg_wc_fst, probs = 0.99, na.rm = TRUE), # Calculate the 99% quantile
         midpoint = (window_pos_1 + 5000)) %>% 
  filter(avg_wc_fst >= quantile_threshold) %>% 
  ungroup() 

#add pop key
fst_peaks <- fst_peaks %>%
  left_join(hr.genome, by = c("chromosome" = "chr")) %>%
  rowwise() %>%
  mutate(fst_pop_key = generate_pop_key(pop1, pop2)) %>%
  ungroup() 

# --- Observed Distances (PST → nearest FST midpoint) ---
observed_distances <- top_1_pct_pst %>%
  inner_join(fst_peaks, by = c("pop_key" = "fst_pop_key", "chr" = "chromosome")) %>%
  mutate(min_distance = abs(midpoint - pos)) %>%
  group_by(pop_key, chr, pos) %>%
  slice_min(min_distance) %>%
  ungroup() %>%
  mutate(distance_normalized = min_distance / max_bp)

# View the top 1% PST data with distances
head(observed_distances)
summary(observed_distances$min_distance)
summary(observed_distances$distance_normalized)

hist(observed_distances$min_distance)
hist(observed_distances$distance_normalized)

#This forms the true dataset with the calculated distances that will be compared to the permutations

# --- Permutation Test ---
#I want to make a random dataset of PST values that randomly pulls X number of permutations, 
#and pulls the same number of PST loci in the same proportions per chromosome and be compared to the FST values 
#to get a distribution of randomized distances.

results.gw <- tibble(
  pop_key = character(),
  observed_mean_distance = numeric(),
  random_mean = numeric(),
  random_std_dev = numeric(),
  random_5 = numeric(),
  random_95 = numeric(),
  mean_compare = logical(),
  p_value = numeric()
)

random_distances_means_gw <- NULL  # DataFrame to store simulation means for all combinations

set.seed(42)  # Set seed for reproducibility

# Unique pop_keys
unique_pop_keys <- unique(observed_distances$pop_key)

# Loop through each pop_key
for (pop_key_vector in unique_pop_keys) {
  # Filter poi_positions for current pop_key
  poi_positions_pop <- observed_distances %>%
    filter(pop_key == pop_key_vector)
  
  popa <- unique(poi_positions_pop$pop1.x)
  popb <- unique(poi_positions_pop$pop2.x)
  two.pop <- unique(poi_positions_pop$pop_key)
  
  #Random datset PST
  row_counts <- observed_distances %>%
    group_by(pop1.x, pop2.x, chr) %>%
    filter(pop1.x == popa, pop2.x == popb) %>%
    summarise(row_count = n(), .groups = "drop")
  
  peak_chr <- fst_peaks %>%
    group_by(pop1, pop2, chromosome) %>%
    filter(fst_pop_key == two.pop) %>%
    distinct(midpoint, max_bp)
  
  peak_chr_list <- peak_chr %>%
    group_by(chromosome) %>%
    summarise(midpoint = list(midpoint), chr_length = list(unique(max_bp)), .groups = "drop")
  
# Filter pst_data_long for random sampling
  random_positions <- pst_data_long %>%
    filter(pop1 == popa, pop2 == popb)
  
  if (nrow(poi_positions_pop) > 0 && nrow(random_positions) > 0) {
    n_pois <- nrow(poi_positions_pop)
    n_simulations <- 1000
    
    random_distances_mean <- numeric(n_simulations)
    
    for (i in 1:n_simulations) {
      # Sample random points, maintaining chromosome distribution
      sample_chr <- function(chr_name) {
        # Subset rows for the given chromosome
        chr_data <- random_positions %>% filter(chr == chr_name)
        # Get the number of rows to sample
        n_samples <- row_counts %>% filter(chr == chr_name) %>% pull(row_count)
        # If no matching count is found, return an empty tibble (this shouldnt happen)
        if (length(n_samples) == 0 || is.na(n_samples)) return(tibble())
        # Perform random sampling, ensuring we don't sample more than available
        chr_data %>% slice_sample(n = min(n_samples, nrow(chr_data)))
      }
      
      # Apply function to all unique chromosomes in random_positions
      sampled_points <- map_dfr(unique(random_positions$chr), sample_chr)
      
      sampled_points <- sampled_points %>%
        left_join(peak_chr_list, by = c("chr" = "chromosome")) %>%
        # Drop chromosomes without midpoint info
        filter(!map_lgl(midpoint, ~ is.null(.x) || length(.x) == 0)) %>%
        mutate(
          distances = map2_dbl(pos, midpoint, ~ {
            chr_midpoints <- unlist(.y)
            min(abs(.x - chr_midpoints), na.rm = TRUE)
          }),
          distances = distances / as.numeric(unlist(chr_length))
        )
      
      # Store the mean distance for this simulation
      random_distances_mean[i] <- mean(sampled_points$distances)
    }

    # Observed distances
    poi_positions_pop <- poi_positions_pop %>%
      mutate(min_distances = min_distance / max_bp)
    observed_mean <- mean(poi_positions_pop$min_distances)
    p_value <- mean(random_distances_mean <= observed_mean)
    
    # Store the main results
    results.gw <- bind_rows(results.gw, tibble(
      pop_key = pop_key_vector,
      random_mean = mean(random_distances_mean),
      random_std_dev = sd(random_distances_mean),
      random_5 = quantile(random_distances_mean, 0.05, na.rm=T),
      random_95 = quantile(random_distances_mean, 0.95, na.rm=T),
      observed_mean_distance = observed_mean,
      mean_compare = observed_mean <= mean(random_distances_mean),
      p_value = p_value
    ))
    
    # Store the random simulation means in a separate dataframe
    temp_means_df <- as.data.frame(t(random_distances_mean))
    colnames(temp_means_df) <- paste0("Sim_", seq_len(n_simulations))
    temp_means_df <- cbind(
      pop_key = pop_key_vector,
      temp_means_df
    )
    random_distances_means_gw <- bind_rows(random_distances_means_gw, temp_means_df)
  }
}

# Join the simulation means with the main results
final_results_gw <- left_join(results.gw, random_distances_means_gw, by = c("pop_key"))

final_results_long <- final_results_gw %>%
  pivot_longer(
    cols = starts_with("Sim_"),  # Pivot all simulation columns
    names_to = "simulation",    # New column for simulation identifiers
    values_to = "simulated_value"  # New column for simulation values
  )


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


write.table(final_results_gw, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_PST_nearest_FST/2_simulation_results/1_PST_99_FST_99_1000sim_all_results.txt", quote=F, row.names=F)
#final_results_gw <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_simulation_results/1_PST_99.9_nearest_FST_99_1000sim_results.txt", header=T)

