#!/bin/bash
#make positions file by using awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1}' meth.168.20miss.txt > meth.168.20miss.positions.bed
#VCF=vcf.gz
METH=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/meth.168.merged.20miss.positions.bed
FILE=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/1_homer/homer.input.168.20miss.bed

#for vcf file
#bcftools query -f '%CHROM\t%POS0\t%END\n' $VCF > new.file.bed
#echo "chr start end id strand" > $FILE
#awk '{print $0, $1"."$3, 1}' new.file.bed >> $FILE
#sed -i 's/ /\t/g' $FILE

#for methylation file (2025)
echo "chr start end id strand" > $FILE
awk '{print $1, $2, $2, $1"_"$2, 1}' $METH >> $FILE
sed -i 's/ /\t/g' $FILE
