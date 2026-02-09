#!/bin/bash

# Input and output file paths
input_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/HR.168.invariant/HR.168.invariant.eu_pi.txt"
output_file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_PI/1_pi_avg_gw.txt"

# Remove header from input file
tail -n +2 "$input_file" | awk '
{
  # Create a key for each unique population
  pop = $1

  # Sum avg_pi values and count occurrences for each population
  pi_sum[pop] += $5
  count[pop] += 1
} 
END {
  # Output the average pi for each population
  for (pop in pi_sum) {
    if (count[pop] > 0) {
      avg_pi = pi_sum[pop] / count[pop]
      print pop "\t" avg_pi
    }
  }
}' > "$output_file"

# Add a header to the output file
echo -e "pop\tavg_pi" | cat - "$output_file" > temp && mv temp "$output_file"

echo "Averaged pi data saved to $output_file"
