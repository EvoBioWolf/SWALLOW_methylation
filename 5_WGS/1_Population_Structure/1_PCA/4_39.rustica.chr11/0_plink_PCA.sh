#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=1:00:00
#SBATCH -J PCA.plink
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH --error=PCA_plink.err
#SBATCH --out=PCA_plink.out

##### Swallow Methylation -39 rustica rustica individuals
#Perform PCA 
#Written by: Sarah Mueller 
#Modification date: 2023-04-20
### Use: sbatch 1_PCA_plink.sh <vcf complete path> <outdir> <output prefix> <prune.in file> <remove ind files>
#___________ PREP _______________

#We need a previously filter vcf file as an input
vcf=$1
outdir=$2
output=$outdir/$3
#prunein=$4
keep=$4
#___________ MAIN _______________

# use LD files prune and create pca
plink --vcf $vcf \
      --keep $keep \
      --chr 11 --from-bp 18200001 --to-bp 21491856 \
      --double-id \
      --allow-extra-chr \
      --allow-no-sex \
      --keep-allele-order \
      --chr-set 37 no-xy \
      --set-missing-var-ids @:# \
      --make-bed \
      --pca \
      --out $output

#Describing options used in the script
<<COMMENT
--extract - conduct the analysis removing prunned sites
--make-bed - text file format used to store genomic regions as coordinates and associated annotations. -- for admixture analysis
--pca - calculate a principal components analysis.
--remove to remove unwanted individuals from analysis
   --chr 11 --from-bp 18200001 --to-bp 21491856 \
     --chr 03 --from-bp 12550001 --to-bp 13800000 \
       --chr 12 --from-bp 8600001 --to-bp 8900000 \
COMMENT
