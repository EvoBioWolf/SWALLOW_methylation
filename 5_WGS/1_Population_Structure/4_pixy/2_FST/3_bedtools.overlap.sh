#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=20:00
#SBATCH -J intersect

#run bedtools to find the CpG sites which overlap with repeat regions and CpG islands
#this will later be combined with homer output to look at genomic architecture of CpG sites
bedtools intersect -wa -wb -nonamecheck \
                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_windows.txt \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged_update.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_cpgi.txt

bedtools intersect -wa -wb -nonamecheck \
                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_windows.txt \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.merged.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_shore.txt

#bedtools intersect -wa -wb -nonamecheck \
#                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_windows.txt \
#                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shelves.bed \
#                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_shelf.txt

bedtools subtract -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_windows.txt \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.merged.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_opensea.tmp

bedtools subtract -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_opensea.tmp \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged_update.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/2_FST/3_fst_opensea.txt

rm *tmp

#Describing options used in the script
<<COMMENT
  -wa keeps the entries from the A input file
  -wb keeps the entries from the B input file, only for overlapping regions
  -nonamecheck For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. "chr1" vs "chr01".  
  -a input file A (the file used to find overlaps in B)
  -b the file used for finding any overlaps from A
COMMENT
