#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=30:00
#SBATCH -J bismark.extraction

#Sarah Mueller --based on script from Justin Merundon
#For filtering merged cov files from methylation extraction in bismark
#Usage: for i in $(ls *M.CpG_merged.cov | rev | cut -c16- | rev); do sbatch -J ${i} /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_bismark.filter.merged.combo_p.sh ${i}; done
#Usage: for i in $(ls *F.CpG_merged.cov | rev | cut -c16- | rev); do sbatch -J ${i} /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_bismark.filter.merged.combo_p.sh ${i}; done
#Usage: for i in $(ls *U.CpG_merged.cov | rev | cut -c16- | rev); do sbatch -J ${i} /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_bismark.filter.merged.combo_p.sh ${i}; done


SCRATCH=/tmp
RUN=$SLURM_JOB_NAME


#Feature files 
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2
cgi=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/HR_CGI_merged.bed 
snps=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/2_CT_filter/all_c_g_snps.bed
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/2_repeatmasker/HR.repeats.bed
chrs=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome
chrsbed=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.bed
sexdmr=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/3_sex_evaluation/sex_dmr_remove.bed
cutsites=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_cut_sites/mspI_5p_chr_seqkit_150bp.bed

#Directory paths 
workdir=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/1_bismark_map/1_output_covfiles
output=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_output_files_test

#Need to copy the input files into the 'SCRATCH' directory (this could be the fastq file, bismark extraction, bam etc.) 
#Since this is for filtering, then its the bismark extraction output that we copy over 
cat ${workdir}/${RUN}.CpG_merged.cov > $SCRATCH/${RUN}.CpG_merged.cov

cd ${SCRATCH}

##Variables
#Minimum Coverage Required for a site to be kept
mincov=5
maxcov=200


#Only keep sites on the major chromosomes 
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}.CpG_merged.cov | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' > ${SCRATCH}/${RUN}.tmp1

#Only keep positions within predicted inserts
bedtools intersect -wb -a ${cutsites} -b ${SCRATCH}/${RUN}.tmp1 | awk '{OFS="\t"}{print $7, $8, $9, $10, $11, $12}' | sort | uniq > ${SCRATCH}/${RUN}.tmp2

#Remove C-T SNP positions on major chromosomes
bedtools subtract -a ${SCRATCH}/${RUN}.tmp2 -b ${snps} > ${SCRATCH}/${RUN}.tmp3

#Remove sex related positions on major chromosomes
bedtools subtract -a ${SCRATCH}/${RUN}.tmp3 -b ${sexdmr} > ${SCRATCH}/${RUN}.tmp4

#Filter for sites with at least 5 reads and no more than 200 reads
awk -v x=${mincov} -v y=${maxcov} '($5 + $6) >= x && ($5 + $6) <= y' ${SCRATCH}/${RUN}.tmp4 > ${SCRATCH}/${RUN}.tmp5

#Filter out those from repeat regions 5mCs
bedtools intersect -wb -a ${repeats} -b ${SCRATCH}/${RUN}.tmp5 | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | gzip -c > ${output}/${RUN}.CpG_in_Repeats.cov.gz
bedtools subtract -a ${SCRATCH}/${RUN}.tmp5 -b ${repeats} | gzip -c > ${output}/${RUN}.CpG_5mC.cov.gz

### SUMMARIZE COUNTS
#Starting, raw positions
echo "RAW" >> ${output}/${RUN}.COUNTS.txt
cat ${workdir}/${RUN}.CpG_merged.cov| wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining on major scaffolds
echo "MAJOR_CHROM" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp1 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after cut site filter
echo "CUT_SITE_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp2 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering SNPs
echo "SNP_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp3 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering sex related positions
echo "SEX_DMR_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp4 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering for coverage
echo "COVERAGE_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp5 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions in repeats
echo "REPEAT_FILTER" >> ${output}/${RUN}.COUNTS.txt
zcat ${output}/${RUN}.CpG_5mC.cov.gz | wc -l >> ${output}/${RUN}.COUNTS.txt
