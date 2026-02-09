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
                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/meth.168.merged.20miss.positions.bed \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged_update.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgisland_cpgloci_overlap/cpgisland_cpgloci_output.168.20miss.update.bed

#Describing options used in the script
<<COMMENT
  -wa keeps the entries from the A input file
  -wb keeps the entries from the B input file, only for overlapping regions
  -nonamecheck For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. "chr1" vs "chr01".  
  -a input file A (the file used to find overlaps in B)
  -b the file used for finding any overlaps from A
COMMENT
