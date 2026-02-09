#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00
#SBATCH -J repeatmasker

GFF=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta.out.gff

awk '{print $1, $4, $5, $7, $10}' $GFF > repeats.tmp1
grep -E -- 'chr' repeats.tmp1 > repeats_chrom.tmp
awk '{OFS="\t"}{print $1, $2, $3, $4, $5}' repeats_chrom.tmp | sed 's/"//g' > repeats.tmp2
cp repeats.tmp2 HR_repeats.bed
