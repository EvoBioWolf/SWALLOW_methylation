#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=10:00:00
#SBATCH -J LDprun.plink
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH --error=LD_prunnning_plink.err
#SBATCH --out=LD_prunnning_plink.out


##### Swallow WGS - RRBS Project
#Prune SNPs by Linkage disequilibrium
#Written by: Sarah Mueller (s.mueller@bio.lmu.de) & Paulina Nunez-Valencia <paulina.nunez.ilp@gmail.com>
#Modification date: 2023-06-19
### Use: sbatch 1_LDprunning_plink.sh <vcf complete path> <outdir> <output prefix> <samples to keep file>

#___________ PREP _______________

#We need a previously filter vcf file as an input
vcf=$1

#As an output we are going to obtain a list of SNPs to be prune
outdir=$2
output=$outdir/$3

#___________ MAIN _______________

# And we perform linkage pruning - i.e. identify prune sites
plink --vcf $vcf \
      --double-id \
      --allow-extra-chr \
      --allow-no-sex \
      --keep-allele-order \
      --chr-set 37 no-xy \
      --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 \
      --out $output 


#Then edit the plink.prune in file to make it a chrom.position file for vcftools
#sed -i 's/:/\t/g' $output.prune.in


#Describing options used in the script
<<COMMENT
--vcf - specified the location of our VCF file.
--double-id causes both family and within-family IDs to be set to the sample ID.
--allow-extra-chr - allow additional chromosomes beyond the human chromosome set. 
--set-missing-var-ids - Whole-exome and whole-genome sequencing results frequently contain variants which have not been assigned standard IDs. If you don't want to throw out all of that data, you'll usually want to assign them chromosome-and-position-based IDs. This would name our variants as: chromosome:position. 
--indep-pairwise - 1 - window of 50 Kb. 
                    2- 10 window step size - meaning we move 10 bp each time we calculate linkage. 
                    3 - we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate. Here we prune any variables that show an r2 of greater than 0.1.
--out Produce the prefix for the output data.
--chr-set 37 no-xy changes the chromosome set.  since there is more chr than in humans 64
--keep-allele-order Use this EVERY SINGLE TIME you call a plink command, otherwise the order of Allele1 and Allele2 may (or probably will) be flipped in your data. \
--allow-no-sex PLINK will default to removing individuals that have unassigned sex, use this to force it to keep them. \

--gzvcf original vcf file
--positions loci that passed linkage pruning and should be kept in the new pruned vcf file
--remove-indels remove sites with indels
--recode used for creating new vcf format, the INFO field is not kept, if you want to keep the INFO field you need to use recode-INFO (see vcftools manual)
--stdout prints the vcf file to the terminal so it can be directly zipped 

COMMENT
