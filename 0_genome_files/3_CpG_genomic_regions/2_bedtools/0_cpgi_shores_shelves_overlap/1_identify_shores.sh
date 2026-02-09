#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=20:00
#SBATCH -J intersect

#identify 200bp around the cpgisland...this will define the shore
bedtools slop -b 2000 -i /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged_update.bed \
                 -g /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/bedtools.cpgi.with.shores.bed

#Subtract all cpgislands, so you are only left with shores
bedtools subtract -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/bedtools.cpgi.with.shores.bed \
                 -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged_update.bed \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.bed

#Sort the cpgi.shores.bed file
sort -k1,1 -k2,2n /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.bed \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.sorted.bed

#Merge any overlapping shores
bedtools merge -i /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.sorted.bed \
               > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.merged.bed

bedtools intersect -wa -wb -nonamecheck \
                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/meth.168.merged.20miss.positions.bed \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.merged.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.168.20miss.update.bed

<<COMMENT
slop
    -b the number of basepairs to keep on either side of the genomic regi-
    -i the file to extend 
    -g the genome file
subtract
    -a file to remove regions from
    -b file of regions to remove
intersect
  -wa keeps the entries from the A input file
  -wb keeps the entries from the B input file, only for overlapping regions
  -nonamecheck For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. "chr1" vs "chr01".  
  -a input file A (the file used to find overlaps in B)
  -b the file used for finding any overlaps from A
COMMENT


