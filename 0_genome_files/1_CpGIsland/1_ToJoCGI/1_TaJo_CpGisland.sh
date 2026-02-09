#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=6:00:00
#SBATCH -J TaJoCGI

python TaJoCGI.py -o /dss/dsshome1/lxc0B/ra52qed/output/0.genome.files/1.ToJoCGI/Hirundo_rustica_CGI.bed /dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta

