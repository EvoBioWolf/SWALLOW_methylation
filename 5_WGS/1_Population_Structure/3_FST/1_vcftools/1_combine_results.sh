#!/bin/bash

# File paths
input_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/3_FST/results/vcftools_fst.err"
output_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/3_FST/results/1_combined_fst_results_test.txt"

# Extract lines starting at line 13
column1=$(awk 'NR >= 13 && NR % 32 == 13' "$input_file")

# Extract lines starting at line 31
column2=$(awk 'NR >= 31 && NR % 32 == 31' "$input_file")

# Extract lines starting at line 32
column3=$(awk 'NR >= 32 && NR % 32 == 0' "$input_file")

# Combine and save the extracted lines in the output file
paste <(echo "$column1") <(echo "$column2") <(echo "$column3") >> "$output_file"

sed -i 's/--out//g' $output_file
sed -i 's/Weir and Cockerham mean Fst estimate://g' $output_file
sed -i 's/Weir and Cockerham weighted Fst estimate://g' $output_file