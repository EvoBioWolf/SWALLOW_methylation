#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6000mb
#SBATCH --time=6:00:00
#SBATCH -J makeCGI

Rscript --vanilla /dss/dsshome1/lxc0B/ra52qed/scripts/0.genome.files/1.makeCGI/makeCGI.R

