#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=5:00
#SBATCH -J simRad.swallows

##21Jul2021 by Sarah Mueller
#rename chromosomes in Hirundo.rustica reference genome
##need a new.chr.txt file that has one name per line corresponding to chromosome order in fasta file
#convert long chromosome names into short
seqtk seq -C Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta > Hir.rus.tmp

#get list of chromosome names
grep "^>" Hir.rus.tmp >> chr.names.txt
#remove >
sed -i 's/>//g' chr.names.txt
#make replacement file for seqkit
paste -d"\t" chr.names.txt /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/new.chr.txt >> replace.txt

#use seqkit to replace all chromsome names
seqkit replace -p "(.+)$" -r "{kv}" -k replace.txt Hir.rus.tmp >> 1Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta

mv 1Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta
rm chr.names.txt
rm replace.txt
rm new.chr.txt

