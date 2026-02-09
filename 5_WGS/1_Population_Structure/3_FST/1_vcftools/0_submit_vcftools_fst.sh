#!/bin/bash

# Path to VCF file
vcf="/dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/HR.168.biallelic.maf05.dp_3_60.vcf.gz"

# List of populations
populations=("R_r" "R_eu" "E" "G" "S" "T" "TV")

# Loop through all unique combinations of populations
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1="${populations[$i]}"
        pop2="${populations[$j]}"
        
        # Define the population file paths
        pop1_file="0_pop_files/${pop1}_pop.txt"
        pop2_file="0_pop_files/${pop2}_pop.txt"
        
        # Define the output prefix
        output_prefix="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/3_FST/1_vcftools/1_results/${pop1}_${pop2}_FST"

        # Submit a job to Slurm
        sbatch 0_weir_fst.sh $vcf $pop1_file $pop2_file $output_prefix
        
        # Wait for 5 seconds before submitting the next job
        sleep 2
    done
done