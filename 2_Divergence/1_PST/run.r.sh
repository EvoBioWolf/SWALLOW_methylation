#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=2:00:00
#SBATCH -J run.r

Rscript -e "rmarkdown::render('/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_PST_all_pop_heatmap.Rmd')"
