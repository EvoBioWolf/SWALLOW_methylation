#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=9:00:00
#SBATCH -J bismark.map

workdir=./trim.galore
genome=./genome/v2
output=./bismark.map.score.min

while read f;
do bismark --parallel 5 --score_min L,0,-0.6 --output_dir ${output} --genome ${genome} -1 ${f}__R1__*val_1* -2 ${f}__R2__*val_2*; 
done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/sample.names.txt
