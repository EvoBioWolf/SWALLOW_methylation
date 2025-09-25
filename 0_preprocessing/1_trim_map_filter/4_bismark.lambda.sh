#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=70
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=36:00:00
#SBATCH -J bismark.map

workdir=./trim.galore
genome=./Escherichia.phage.Lambda
output=./mapping.efficency


#bismark_genome_preparation ${genome}

while read f;
do bismark --parallel 5  --output_dir ${output} --genome ${genome} -1 ${f}__R1__*val_1* -2 ${f}__R2__*val_2*; 
echo ${f}; awk 'NR==9 {print}' ${f}*report*;
done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/sample.names.full.txt

cd ${output}
bismark_methylation_extractor -p --no_overlap --bedgraph --cytosine_report --output ${output}/methylation.extraction --parallel 8 --buffer_size 20G --genome_folder ${genome} *bismark_bt2_pe.bam
