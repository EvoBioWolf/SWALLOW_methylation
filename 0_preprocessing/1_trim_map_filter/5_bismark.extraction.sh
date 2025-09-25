#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=6000mb
#SBATCH --time=16:00:00
#SBATCH -J bismark.extraction

#SCRATCH=/tmp
#RUN=$SLURM_JOB_NAME

#p is for paired end
##no_overlap is to only read CpG sites from overlapping R1 and R2 once

bismark_methylation_extractor -p --no_overlap --bedgraph --cytosine_report --output ./bismark.map.score.min --multicore 6 --buffer_size 20G --genome_folder ./genome/v2/ DE_08_R_ZD007739_BL_YRL_M.merged.bam

