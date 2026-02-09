#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=23:00:00
#SBATCH -J run.r

Rscript --vanilla /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/2_Genome_Scan/2_PST_nearest_FST/1_PST_near_FST_simulations_cpgi_window.r