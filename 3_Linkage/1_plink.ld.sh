#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=1:00:00
#SBATCH -J LDprun.plink
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH --error=LD_prunnning_plink.err
#SBATCH --out=LD_prunnning_plink.out


##### Swallow WGS - RRBS Project
#Linkage disequilibrium evaluation using PLINK
#Written by: Sarah Mueller (s.mueller@bio.lmu.de) 
#Modification date: 2025.feb
### Use: sbatch 1_LDprunning_plink.sh <prefix to ped and map file> <outdir> <output prefix> 

#___________ PREP _______________

#We the joined methylation and snp .ped and .map file prefix
myfile=$1

#As an output we are going to obtain a list of SNPs to be prune
outdir=$2
output=$outdir/$3
keep=$4

#___________ MAIN _______________

# And we perform linkage pruning - i.e. identify prune sites
plink --ped $myfile.ped \
      --map $myfile.map \
      --allow-extra-chr \
      --allow-no-sex \
      --keep-allele-order \
     --maf 0.2 \
      --chr 11 \
      --chr-set 37 no-xy \
      -r2 gz --ld-window 100 --ld-window-kb 70 \
      --ld-window-r2 0 \
      --make-bed --out $output 


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
--thin - thin randomly thins out the data - i.e. it randomly retains p proportion of the data. Here we set that to 0.1 or 10%. This means that each site has a 10% probability of remaining in the dataset. This is done to ensure the analysis runs quickly and our output isnt too big as it can quickly get out of hand!
--r2 - finally were on to the options for LD! This tells plink to produce squared correlation coefficients. We also provide the argument gz in order to ensure the output is compressed. This is very important as it is easy to produce EXTREMELY large files.
--ld-window - In this case, we set it to 100 bp - i.e. any sites with > 100 sites between them are ignored.
--ld-window-kb - this is the upper end of the LD window. Here we set it to 1000, meaning that we ignore any two sites more than 1 Mb apart in the genome.
--ld-window-r2 - the final LD command - this sets a filter on the final output but we want all values of LD to be written out, so we set it to 0.
--make-bed - this just makes a plink bedfile for future analyses
--out Produce the prefix for the output data.

when running chromosome 11
--chr 11 --from-bp 18200001 --to-bp 21491856 \ ##add this line is needed to run locally on inversion region...also window size must be increased
-r2 gz --ld-window 100000 --ld-window-kb 3500 \
#PARAMETERS CHOSEN FROM PARAMETER TESTING ON SNP DATA FROM WGS

COMMENT


