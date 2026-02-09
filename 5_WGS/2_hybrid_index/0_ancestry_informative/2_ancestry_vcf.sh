#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=30:00
#SBATCH -J filter

#vctools to filter SNP sites from merged vcf file to only keep ancestry informative SNPs
#see https://vcftools.sourceforge.net/man_latest.html for more info on commands
#--max-missing, defined between 0 and 1, 0 is completely missing sites allowed, 1 is no missing sites allowed. Because we have 2 different 'batches' (Safran vcf and Jan vcf), I will set to 90% so that sites that are completely absent in Jans samples do not make it to the final vcf file
#--exclude-bed was used to remove regions in reptitive regions of the genome from /dss/dsshome1/lxc0B/ra52qed/scripts/0.genome.files/2.repeatmasker/repeat.masker.sh file which was made using repeatmasker


vcftools --gzvcf /dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/HR.168.biallelic.maf05.dp_3_60.vcf.gz \
        --positions $1 \
        --recode --stdout | bgzip -c > $2

