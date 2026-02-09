#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=250Gb
#SBATCH --time=72:00:00
#SBATCH -J run.r

Rscript --vanilla /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/3_sex_evaluation/1_DSS_sex_obj.r
