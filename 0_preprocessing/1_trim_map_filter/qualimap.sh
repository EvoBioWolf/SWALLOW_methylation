#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=60:00:00
#SBATCH -J qualimap

while read f;

	do qualimap bamqc --java-mem-size=4G -nt 8 -bam ${f}*.sorted.bam -outdir ./bismark.1102/${f};


done < sample.names.txt
