#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5GB
#SBATCH --time=4:00:00
#SBATCH -J fst_vcftools


#Script to calculate FST from vcf file
#Usage:need text files that have all samples from each population (one individual per line) called [pop name]_pop.txt and define the populations in the bash script 1_submit_vcftools_fst.sh
#Run 1_submit_vcftools_fst.sh to submit jobs to slurm

#_________________PREP_________________
vcf=$1
pop1=$2
pop2=$3
output_prefix=$4

#________________MAIN__________________
# Run vcftools for the given population pair
vcftools --gzvcf $vcf --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out $output_prefix