#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=18:00:00
#SBATCH -J trim.galore

#-q removes low quality bases (<20) first for RRBS samples, then moves onto adapter removal
#--fastqc runs fastqc on the resulting .fastq file
#--gzip compresses output files
#--retain_unpaired keeps unpaired reads that one read was removed due to being too short after adapter removal/trimming
#--rrbs and --paired for our specific data

while read f;
	do trim_galore --rrbs --paired --clip_r1 8 --clip_r2 10 --quality 25 --fastqc --fastqc_args "--outdir /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2024_SwallowTechreps/99_test_trimgalore" --gzip --retain_unpaired --output_dir /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2024_SwallowTechreps/99_test_trimgalore --cores 4 ${f}__R1__* ${f}__R2__*;
	done < /dss/dsshome1/lxc0B/ra52qed/sample.names.test.txt


