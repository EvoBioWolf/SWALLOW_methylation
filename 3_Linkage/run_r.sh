#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=50GB
#SBATCH --time=2:00:00
#SBATCH -J run.r

Rscript --vanilla /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_convert_meth_to_ped.R