#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=5G
#SBATCH --time=18:00:00
#SBATCH -J newhybrids_parallel

Rscript --vanilla /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RT_hyb/1_newhybrids_script.R