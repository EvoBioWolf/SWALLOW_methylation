#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=60:00
#SBATCH -J vcf2bed

#Extract C-T and G-A and subtract one position from the G-A SNPs since we merge before filtering so these will show up as bp-1
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'REF="C" | REF="G"' /dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/HR.168.filtered.order.vcf.gz | awk '{
  if ($3 == "G") {
    print $1 "\t" ($2 - 1)
  } else if ($3 == "C") {
    print $1 "\t" $2
  }
}' > HR.168.filtered.allsnp.tmp


#output all positions from vcf file
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'REF="C" | REF="G"' /dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/Safran_lab/hirundo.rustica.300.id.allsites.renamed.miss02.maf05.ingroup.vcf.gz | awk '{
  if ($3 == "G") {
    print $1 "\t" ($2 - 1)
  } else if ($3 == "C") {
    print $1 "\t" $2
  }
}' > HR.300.m02.m05.d2_30.allsnps.tmp

#Add the third column to the text file to create a bed file
awk '{print $1 "\t" $2 "\t" $2}' HR.168.filtered.allsnp.tmp > HR.168.filtered.allsnp.bed

#take first 2 columns and convert to 0-based format
awk '{print $1 "\t" $2 "\t" $2}' HR.300.m02.m05.d2_30.allsnps.tmp > HR.300.m02.m05.d2_30.allsnp.bed

#Concatenate the two files to create one bed file for extraction
cat HR.168.filtered.allsnp.bed HR.300.m02.m05.d2_30.allsnp.bed > all_c_g_snps.bed
