#!/bin/bash

#RUN IN COMMAND LINE - For reference sequence
#bedtools makewindows -g /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome -w 10000 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/5_GC_content/1_windows.bed

#bedtools nuc -fi /dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta -bed  /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/5_GC_content/1_windows.bed > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/5_GC_content/1_gc_content.txt

#RUN FOR SUBSPECIES GC CONTENT
# Directory containing subspecies genomes
genome_dir="/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/subspecies_specific"

# BED file for windows
windows_bed="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/5_GC_content/1_windows_5kb.bed"

# Output directory for GC content results
output_dir="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/5_GC_content"
mkdir -p "$output_dir"

# Loop through each subspecies genome (FASTA files)
for genome_fasta in "$genome_dir"/*.fa; do
    # Extract subspecies name from filename
    subspecies=$(basename "$genome_fasta" .fa)
    
    # Define output file for GC content
    output_file="$output_dir/${subspecies}_gc_content.txt"
    
    # Calculate GC content using bedtools nuc
    echo "Calculating GC content for $subspecies..."
    bedtools nuc -fi "$genome_fasta" -bed "$windows_bed" > "$output_file"
    
    echo "GC content saved to $output_file"
done

echo "GC content calculation completed for all subspecies."
