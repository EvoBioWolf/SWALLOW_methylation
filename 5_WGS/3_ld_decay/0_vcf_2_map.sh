#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --time=4:00:00
#SBATCH -J vcf.2.map
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb

##### Swallow WGS - RRBS Project
#Prune SNPs by Linkage disequilibrium
#Written by: Sarah Mueller (s.mueller@bio.lmu.de)

#Create plink file
#plink --file /dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/HR.168.biallelic.maf05.dp_3_60.vcf.gz --allow-extra-chr --recode --out 0_HR_168_maf05_dp3_60

#The plink .ped and .map files from vcftools does not match the exact format we need. We need to remove empty space to have the two alleles from each loci sit together
#awk '{
#    printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6;
#    for (i = 7; i <= NF; i += 2) {
#        printf "\t%s%s", $i, $(i+1);
#    }
#    printf "\n";
#}' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/3_ld_decay/0_HR_168_plink.ped > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/3_ld_decay/0_HR_168_plink.format.ped

#Fix characters in sample names in methylation ped
sed -i 's/2171./2171-/g' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.ped
sed -i 's/3121./3121-/g' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.ped

#First append every line from the methylation PED file to the SNP PED file based on sample name (in my case they are not in the same sample order)
awk 'NR==FNR {for (i=6; i<=NF; i++) data[$1] = data[$1] FS $i; next} $1 in data {print $0, data[$1]}' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.ped /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/3_ld_decay/0_HR_168_plink.format.ped > 0_HR.85.20miss.rmlowv.rmout.snp.meth.ped

# You can double check via the following command:
#awk '{print NF; exit}' # this will print the number of fields in each line, make sure new ped file is equal to sum of old ped files minus the intro fields
#head -n 1 0_HR_168_plink.ped | cut -c1-10 # use this to look at the start of the second file and grep for the name of individuals to check a few and make sure they were assigned correctly...you need to strat at the field where the first file ends

#Next re-configure the MAP file to contain the same chromosome names and identifiers for SNP (S_) or methylation (M_) loci
#First SNP map file
#awk '{
#    split($2, arr, ":");  # Split column 2 at ":"
#    $1 = arr[1];          # Replace column 1 with extracted chromosome part
#    sub(/^chr/, "S_chr", $2);  # Add "S_" before "chr" in column 2
#    print
#}'  /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/3_ld_decay/0_HR_168_plink.map > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/3_ld_decay/0_HR_168_plink.map.txt

#Second methylation map file
awk '{ $2 = "M_" $2; print }' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.map > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.temp.map

#Now merge the two map files
cat /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/3_ld_decay/0_HR_168_plink.map /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.temp.map > 0_HR.85.20miss.rmlowv.rmout.snp.meth.map

#Remove temporary file 
rm  /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.temp.map
