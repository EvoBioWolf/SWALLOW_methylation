#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=10:00
#SBATCH -J bismark.map

while read f;

	do echo ${f}; awk 'NR==9 {print}' ${f}*report*;
	done < sample.names.full.txt
