# Load the necessary library
library(readr)
library(dplyr)

# Read the file (no header), specifying the column names
cpgi_merged <- read_delim(
  "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged.bed",
  delim = "\t", # Assuming the file is tab-separated. Use " " if space-separated.
  col_names = c("chr", "start", "end"),
  show_col_types = FALSE
)

# Create the new 'length' column by subtracting 'start' from 'end'
cpgi_merged <- cpgi_merged %>%
  mutate(length = end - start)

# summarizes the lengths of cpgi islands
summary(cpgi_merged$length)


# Read the file (no header), specifying the column names
cpgi_tajo <- read_delim(
  "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/1_ToJoCGI/Hirundo_rustica_TaJO_CGI.bed",
  delim = "\t", # Assuming the file is tab-separated. Use " " if space-separated.
  col_names = c("chr", "start", "end"),
  show_col_types = FALSE
)

# Create the new 'length' column by subtracting 'start' from 'end'
cpgi_tajo <- cpgi_tajo %>%
  mutate(length = end - start)

# summarizes the lengths of cpgi islands
summary(cpgi_tajo$length)

# Read the file (no header), specifying the column names
cpgi_make <- read_delim(
  "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/1_makeCGI/Hirundo_rustica_makecgi_250.bed",
  delim = "\t", # Assuming the file is tab-separated. Use " " if space-separated.
  col_names = c("chr", "start", "end"),
  show_col_types = FALSE
)

# Create the new 'length' column by subtracting 'start' from 'end'
cpgi_make <- cpgi_make %>%
  mutate(length = end - start)

# summarizes the lengths of cpgi islands
summary(cpgi_make$length)

###############
##I realized that CpG islands can be very short because i used bedtools intersect where it only keeps the overlapping portion, when I should have kept the whole CpG island that was called
##Now I went back and used the -wao option to merge them differently using bedtools intersect

cpgi_merge_long <- read_delim(
  "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged_update.bed",
  delim = "\t", # Assuming the file is tab-separated. Use " " if space-separated.
  col_names = c("chr", "start", "end"),
  show_col_types = FALSE
)

cpgi_merge_long <- cpgi_merge_long %>% filter(!X4 == ".")
cpgi_merge_long <- cpgi_merge_long %>% mutate(distance = end - start)


#calculate distances between cpgs (i get some overlapping shores, which I want to avoid)
cpgi_merge_long <- cpgi_merge_long %>%
  mutate(distance_to_next = lead(start) - end)
