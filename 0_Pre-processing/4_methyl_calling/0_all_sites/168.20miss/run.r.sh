#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=4:00:00
#SBATCH -J run.r

Rscript --vanilla /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/1_methylkit.unite.20miss.r

#Rscript -e "rmarkdown::render('/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/2_count.20miss.Rmd')"

#Rscript -e "rmarkdown::render('/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/2_rm_outliers_invariant.Rmd')"
