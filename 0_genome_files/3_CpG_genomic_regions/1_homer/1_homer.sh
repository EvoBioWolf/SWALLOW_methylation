#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00
#SBATCH -J homer

##### Swallow WGS - RRBS Project
#Perform homer analysis to get genomic annotation of CpG sites within the HR genome
#Written by: Sarah Mueller <s.mueller@bio.lmu.de>
#Modification date: 2023-08-04
### Use: sbatch 1_homer.sh <ref genome> <gtf file> <input bed file> <output txt file> 

#NOTE: you need to make the CpG loci that you want to be annotated into a bed file for input into homer...most errors arisign are likely due to format of this input file

#___________ PREP _______________

NAME=bHirRus1.pri.v2
ORG=Hirundo_rustica
REF=$1
GTF=$2
INPUT=$3
OUTPUT=$4

#____________MAIN_________________

loadGenome.pl -name $NAME -org $ORG -fasta $REF -gtf $GTF

annotatePeaks.pl $INPUT $REF -gtf $GTF > $OUTPUT


#Describing options used in the script
<<COMMENT
This works best with the gtf file...gff files also supported, but kept giving errors when running with gff

-keepAll keeps all peak sites even if no annotation is found


COMMENT