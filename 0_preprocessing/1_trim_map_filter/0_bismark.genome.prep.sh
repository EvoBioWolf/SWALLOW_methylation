#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00
#SBATCH -J bismark.genome

bismark_genome_preparation /dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2
