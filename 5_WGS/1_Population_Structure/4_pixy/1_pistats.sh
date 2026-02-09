#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=31
#SBATCH --time=2-00:00:00
#SBATCH -J pistats

##### Swallow WGS 
#Calculate nucleotide statistics
#Written by: Sarah Mueller
#Modification date: 2023-01-07
### Use: sbatch 2_pistats_pixy.sh <workdir> <vcf in> <popmap> <out file> 

#___________ PREP _______________

workdir=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/HR.168.invariant.10kb
vcfin=/dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/invariant/3_HR.168.invariant.filtered.order.vcf.gz
popmap=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/1_168_pop_map_eu.txt
out=HR.168.invariant.10kb

#___________ MAIN _______________

#source activate pixy

cd $workdir

# Calculate stats using pixy
pixy --stats pi fst dxy \
     --vcf $vcfin \
     --n_cores 30 \
     --populations $popmap \
     --window_size 10000 \
     --bypass_invariant_check 'yes' \
     --output_folder /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/1_Population_Structure/4_pixy/HR.168.invariant.10kb \
     --output_prefix $out
