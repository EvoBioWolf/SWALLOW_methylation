#!/bin/bash

# Input and output file paths
input_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/HR.168.invariant.5kb/HR.168.invariant.5kb_dxy.txt"
output_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_DXY/1_dxy_avg_5kb_gw.txt"

# Remove header from input file
tail -n +2 "$input_file" | awk '
{pop_pair = $1 "\t" $2

  # Initialize values in associative arrays
  dxy_sum[pop_pair] += $6
  count[pop_pair] += 1
} 
END {
  # Output the averages for each population pair
  for (pair in dxy_sum) {
    if (count[pair] > 0) {
      avg_dxy = dxy_sum[pair] / count[pair]
      print pair "\t" avg_dxy
    }
  }
}' > "$output_file"

# Add a header to the output file
echo -e "pop1\tpop2\tavg_dxy" | cat - "$output_file" > temp && mv temp "$output_file"

echo "Averaged data saved to $output_file"
