##Calculate sequence counts per samples for Muelelr et al. 2025 - Methylation

#From mutliqc data, read in number of unique reads per sample
sequence_counts <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/0_quality/0_rawdata/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt", header=T, sep= "\t")

sequence_counts <- sequence_counts %>%
  separate(
    col = Sample,
    into = c("sample_name", "run"),
    sep = "__S"
 )

#Add together different runs of each sample to get total values
sequence_counts_ind <- sequence_counts %>%
  group_by(sample_name) %>%
  summarise(
    Total_Unique_Reads = sum(`Unique.Reads`),
    Total_Duplicate_Reads = sum(`Duplicate.Reads`)
  ) %>%
  ungroup() # It's good practice to ungroup after summarizing

# Print the final summary data frame
cat("\nAverage Reads per Library:\n")
summary(sequence_counts_ind)
