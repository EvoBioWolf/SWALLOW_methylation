# Install and load necessary packages if you haven't already
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("IRanges", quietly = TRUE)) {
  install.packages("IRanges")
}
library(tidyverse)
library(IRanges) # For efficient interval operations

# 1. Load the data
# Adjust 'header = TRUE' if your files have headers, which they do based on your example
# Adjust 'sep = " "' if your columns are space-separated, or 'sep = "\t"' if tab-separated
dxy_data <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/3_high_dxy_regions.txt", header = TRUE, sep = " ")
fst_data <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/3_high_fst_regions.txt", header = TRUE, sep = " ")

# 2. Prepare the data for comparison
# Ensure column names are consistent and suitable for IRanges
# As per our previous conversation, I'll use dplyr::rename here.
dxy_data <- dxy_data %>%
  dplyr::rename(
    start = window_pos_1,
    end = window_pos_2
  )

fst_data <- fst_data %>%
  dplyr::rename(
    start = window_pos_1,
    end = window_pos_2
  )

# 3. Find overlapping windows for the same pop_key
# This involves iterating through each pop_key and chromosome
overlapping_windows <- list()

# Get unique population keys from both datasets
all_pop_keys <- unique(c(dxy_data$pop_key, fst_data$pop_key))

for (pop in all_pop_keys) {
  # Filter data for the current pop_key
  dxy_pop <- dxy_data %>% filter(pop_key == pop)
  fst_pop <- fst_data %>% filter(pop_key == pop)
  
  # Get unique chromosomes for the current pop_key
  chromosomes <- unique(c(dxy_pop$chromosome, fst_pop$chromosome))
  
  for (chr in chromosomes) {
    # Filter data for the current chromosome
    dxy_pop_chr <- dxy_pop %>% filter(chromosome == chr)
    fst_pop_chr <- fst_pop %>% filter(chromosome == chr)
    
    # Check if both datasets have data for this pop/chr combination
    if (nrow(dxy_pop_chr) > 0 && nrow(fst_pop_chr) > 0) {
      # Create IRanges objects for efficient overlap finding
      dxy_ranges <- IRanges(start = dxy_pop_chr$start, end = dxy_pop_chr$end)
      fst_ranges <- IRanges(start = fst_pop_chr$start, end = fst_pop_chr$end)
      
      # Find overlaps
      overlaps <- findOverlaps(dxy_ranges, fst_ranges)
      
      if (length(overlaps) > 0) {
        # Extract the overlapping windows from both original data frames
        # and combine them. We'll include both original rows for context.
        dxy_overlapping_rows <- dxy_pop_chr[queryHits(overlaps), ]
        fst_overlapping_rows <- fst_pop_chr[subjectHits(overlaps), ]
        
        # You might want to combine these into a single data frame
        # For simplicity, let's create a data frame showing the overlap details
        # You could also join the original rows if you need all columns.
        overlap_details <- data.frame(
          chromosome = chr,
          pop_key = pop,
          dxy_window_start = dxy_overlapping_rows$start,
          dxy_window_end = dxy_overlapping_rows$end,
          fst_window_start = fst_overlapping_rows$start,
          fst_window_end = fst_overlapping_rows$end,
          stringsAsFactors = FALSE
        )
        
        # Calculate the actual overlap coordinates (intersection)
        # This will give you the precise genomic region of overlap
        overlap_ranges <- pintersect(dxy_ranges[queryHits(overlaps)], fst_ranges[subjectHits(overlaps)])
        overlap_details$overlap_start <- start(overlap_ranges)
        overlap_details$overlap_end <- end(overlap_ranges)
        overlap_details$overlap_width <- width(overlap_ranges)
        
        overlapping_windows[[paste(pop, chr, sep = "_")]] <- overlap_details
      }
    }
  }
}

# Combine all results into a single data frame
if (length(overlapping_windows) > 0) {
  final_overlaps_df <- do.call(rbind, overlapping_windows)
  rownames(final_overlaps_df) <- NULL # Clean up row names
  print("Overlapping windows found:")
  print(final_overlaps_df)
} else {
  print("No overlapping windows found for any pop_key.")
}

# You can now further analyze 'final_overlaps_df'
# For example, to save the results:
# write.csv(final_overlaps_df, "overlapping_windows_results.csv", row.names = FALSE)