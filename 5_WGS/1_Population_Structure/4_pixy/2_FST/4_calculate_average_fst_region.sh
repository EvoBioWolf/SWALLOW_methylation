#!/bin/bash

# Input and output file paths
input_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_opensea.txt"
output_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/4_fst_opensea_pop.txt"

# Remove header from input file
tail -n +2 "$input_file" | awk '
{pop_pair = $4 "\t" $5

  # Initialize values in associative arrays
  fst_sum[pop_pair] += $6
  count[pop_pair] += 1
} 
END {
  # Output the averages for each population pair
  for (pair in fst_sum) {
    if (count[pair] > 0) {
      avg_fst = fst_sum[pair] / count[pair]
      print pair "\t" avg_fst
    }
  }
}' > "$output_file"

# Add a header to the output file
echo -e "pop1\tpop2\tavg_fst" | cat - "$output_file" > temp && mv temp "$output_file"

echo "Averaged data saved to $output_file"



##Combine all into one output
# paste 4_fst_cpgi_pop.txt <(cut -f 3 4_fst_shore_pop.txt) <(cut -f 3 4_fst_opensea_pop.txt) > 4_fst_pop_combined.txt