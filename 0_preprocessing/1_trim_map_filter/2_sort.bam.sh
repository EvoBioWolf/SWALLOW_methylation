#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=10:00:00
#SBATCH -J sort.bam

while read f;

	do samtools sort -@ 8 -o ${f}.sorted.bam ${f}*.bam;

done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/sample.names.txt
