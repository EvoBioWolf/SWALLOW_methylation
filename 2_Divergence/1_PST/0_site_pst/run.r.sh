#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-12:00:00
#SBATCH -J pstat

Rscript --vanilla /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_site_pst/0_PST_allsite_variant.r