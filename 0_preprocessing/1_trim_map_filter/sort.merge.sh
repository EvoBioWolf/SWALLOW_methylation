#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00
#SBATCH -J sort.merge

SCRATCH=/tmp
RUN=$SLURM_JOB_NAME

#Feature files 
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2

#Directory paths 
workdir=./bismark.map.score.min
output=./bismark.filter

cd ${workdir}

#sort bam file
for i in $(ls *_val_1_bismark_bt2_pe.bam | rev | cut -c26- | rev); 
	do samtools sort ${i}_val_1_bismark_bt2_pe.bam -o ${i}.sorted.bam;
	done


