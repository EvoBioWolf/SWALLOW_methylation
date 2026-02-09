# Load necessary libraries
library(dplyr)
library(readr) # For faster reading of large files

# Define the window size in base pairs (10kb = 10,000 bp)
window_size_bp <- 10000

# --- Read Data ---
# Read the data from the text file.
# For space-separated: read_delim(input_file_path, delim = " ", trim_ws = TRUE)
# For tab-separated: read_tsv(input_file_path)
# For comma-separated: read_csv(input_file_path)
# Using read_delim with whitespace as delimiter, and trimming extra whitespace.
data <- read.table(gzfile("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/HR_85_snp.meth.maf2.ldwin100.gz"), sep="", header=T)
colnames(data) <- c("CHR_A",	"BP_A",	"SNP_A",	"CHR_B",	"BP_B",	"SNP_B",	"R2")

# Check the dataframe to make sure its the meth-meth, snp-meth or snp-snp file
all_snp_a_start_with_M <- all(startsWith(data$SNP_A, "S_"))
print(paste("Do all SNP_A values start with 'M_'?", all_snp_a_start_with_M))
# Check if all SNP_B values start with "M_"
all_snp_b_start_with_M <- all(startsWith(data$SNP_B, "M_"))
print(paste("Do all SNP_B values start with 'M_'?", all_snp_b_start_with_M))

# Ensure BP_A and R2 are numeric
data$BP_A <- as.numeric(data$BP_A)
data$R2 <- as.numeric(data$R2)
data$BP_B <- as.numeric(data$BP_B)
data$CHR_A <- as.numeric(data$CHR_A)

# --- Process Data into Windows ---
# Create 10kb windows and calculate the average R2
averaged_r2_windows <- data %>%
  # Group by chromosome (CHR_A)
  dplyr::group_by(CHR_A) %>%
  # Calculate the start of the 10kb window for each BP_A
  # floor((BP_A - 1) / window_size_bp) * window_size_bp + 1
  # This formula ensures that 1-10000 falls into window 1, 10001-20000 into window 2, etc.
  dplyr::mutate(
    window_start_bp = floor((BP_A - 1) / window_size_bp) * window_size_bp + 1,
    window_end_bp = window_start_bp + window_size_bp - 1
  ) %>%
  # Group by chromosome and the newly created window start
  dplyr::group_by(CHR_A, window_start_bp, window_end_bp) %>%
  # Calculate the average R2 for each window
  dplyr::summarise(
    average_R2 = mean(R2, na.rm = TRUE), # na.rm = TRUE to ignore NA values in R2
    .groups = "drop" # Drop the grouping after summarising
  ) %>%
  # Arrange the results for better readability
  dplyr::arrange(CHR_A, window_start_bp)

# --- Display Results ---
# Print the first few rows of the result
print("Calculated average R2 in 10kb windows:")
print(head(averaged_r2_windows))

# --- Save Results to File ---
# Save the full averaged_r2_windows data to a text file
# 'sep = "\t"' for tab-separated, 'row.names = FALSE' to prevent writing row numbers.
write.table(averaged_r2_windows, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/6_snp_meth_r2_10kb_windows.txt", sep = "\t", row.names = FALSE, quote = FALSE)
print("Full results saved to '6_snp_meth_r2_10kb_windows.txt'")

# --- Filter for chr11 and Save ---
# Filter the averaged_r2_windows for only entries on chromosome "chr11"
chr11_averaged_r2_windows <- averaged_r2_windows %>%
  dplyr::filter(CHR_A == "chr11")
print(head(chr11_averaged_r2_windows))
# Save the chr11-specific results to a separate text file
write.table(averaged_r2_windows, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/6_snp_r2_10kb_windows_chr11.txt", sep = "\t", row.names = FALSE, quote = FALSE)
print("Chr11 results saved to '6_averaged_r2_10kb_windows_chr11.txt'")

