#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=5:00
#SBATCH -J mergeCGI


#taxa=Hirundo_rustica
#Output raw makecgi File
#makecgi=/dss/dsshome1/lxc0B/ra52qed/scripts/0.genome.files/1.makeCGI/CGI-Hirundo_rustica.txt
#Minimum cpg island size 
#makecgiwin=250

#Filter for windows > window size 
#awk -v w=${makecgiwin} '$4 > w' ${makecgi} > ${taxa}_makecgi_${makecgiwin}.tmp

#Output first three columns (chr, start, end) 
#awk '{OFS="\t"}{print $1, $2, $3}' ${taxa}_makecgi_${makecgiwin}.tmp > ${taxa}_makecgi_${makecgiwin}.bed

#remove temporary file
#rm *.tmp

#delete first line in makecgi .bed because it only contains identifiers
#sed -i '1d' ${taxa}_makecgi_${makecgiwin}.bed

#Find the consistent islands between tajo's algorithm and makecgi
#bedtools intersect -woa ${taxa}_makecgi_${makecgiwin}.bed -b ${taxa}_CGI.bed > ${taxa}_CGI_merged.bed

#makeCGI island count
#wc -l Hirundo_rustica_makecgi_250.bed

#Tajo CpG island count
#wc -l Hirundo_rustica_CGI.bed

#Intersected results
#wc -l Hirundo_rustica_CGI_merged.bed

#Find the consistent islands between tajo's algorithm and makecgi, keeping the full island if overlap is found
#bedtools intersect -wao -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/1_makeCGI/Hirundo_rustica_makecgi_250.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/1_ToJoCGI/Hirundo_rustica_TaJO_CGI.bed > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/Hirundo_rustica_CGI_full_merged.bed
#############REMOVE INTERSECT ABOVE IF YOU CHANGE TO UPDATED VERSION BELOW

####UPDATE 2025
##intersect leaves me wiht some CpG islands that are very small. I end up with 280 CpG islands with less than 200bp. These are not very biologically meaningful. I thought a better way than looking for overlaps
#was to merge the two so that islands are consolidated if there are very small gaps between them or overlaps. 

cat /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/1_makeCGI/Hirundo_rustica_makecgi_250.bed /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/1_ToJoCGI/Hirundo_rustica_TaJO_CGI.bed | bedtools merge -d 100 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged_update.bed