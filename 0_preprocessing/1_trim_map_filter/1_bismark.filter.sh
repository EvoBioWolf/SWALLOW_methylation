#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=30:00
#SBATCH -J bismark.extraction

SCRATCH=/tmp
RUN=$SLURM_JOB_NAME


#Feature files 
genome=./genome/v2
cgi=HR_CGI_merged.bed
snps=HR.300.m02.m05.d2_30.bed
repeats=HR.repeats.bed
chrs=chromosomes.genome
chrsbed=chromosomes.bed
cutsites=mspI_5p_chr_seqkit_150bp.bed

#Directory paths 
workdir=./1_bismark_map
output=./2_bismark_filter/1_merged_CpG_covfiles

#Need to copy the input files into the 'SCRATCH' directory (this could be the fastq file, bismark extraction, bam etc.) 
#Since this is for filtering, then its the bismark extraction output that we copy over 
cat ${workdir}/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov.gz > $SCRATCH/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov.gz

cd ${SCRATCH}

##Variables
#Minimum Coverage Required for a site to be kept
mincov=5
maxcov=200


#Only keep sites on the major chromosomes 
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov.gz | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' > ${SCRATCH}/${RUN}.tmp1

#Only keep positions within predicted inserts
bedtools intersect -wb -a ${cutsites} -b ${SCRATCH}/${RUN}.tmp1 | awk '{OFS="\t"}{print $7, $8, $9, $10, $11, $12}' | sort | uniq > ${SCRATCH}/${RUN}.tmp2

#Remove C-T SNP positions on major chromosomes
bedtools subtract -a ${SCRATCH}/${RUN}.tmp2 -b ${snps} > ${SCRATCH}/${RUN}.tmp3

#Filter for sites with at least 5 reads and no more than 200 reads
awk -v x=${mincov} -v y=${maxcov} '($5 + $6) >= x && ($5 + $6) <= y' ${SCRATCH}/${RUN}.tmp3 > ${SCRATCH}/${RUN}.tmp4

#Filter out those from repeat regions 5mCs
bedtools intersect -wb -a ${repeats} -b ${SCRATCH}/${RUN}.tmp4 | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | gzip -c > ${output}/${RUN}.CpG_in_Repeats.cov.gz
bedtools subtract -a ${SCRATCH}/${RUN}.tmp4 -b ${repeats} | gzip -c > ${output}/${RUN}.CpG_5mC.cov.gz

### SUMMARIZE COUNTS
#Starting, raw positions
echo "RAW" >> ${output}/${RUN}.COUNTS.txt
zcat ${workdir}/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov.gz | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining on major scaffolds
echo "MAJOR_CHROM" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp1 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after cut site filter
echo "CUT_SITE_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp2 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering SNPs
echo "SNP_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp3 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering for coverage
echo "COVERAGE_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp4 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions in repeats
echo "REPEAT_FILTER" >> ${output}/${RUN}.COUNTS.txt
zcat ${output}/${RUN}.CpG_5mC.cov.gz | wc -l >> ${output}/${RUN}.COUNTS.txt
