## Ancestry Informative SNPS

These scripts are quite simple. I want to use introgress to define hybrid indices based on ancestry informative SNPs (FST > 0.8). In order to do that I first want to identify high FST SNPs, and then pull them into a vcf file (unpacking the whole vcf file and calling it in R takes much longer...usually I don't like having suplicated information, but here it really helps with speed.)

1_ancestry_informative.sh : finds sites from vcftools FST output (which was already calculated) that are above 0.8

2_ancestry_vcf.sh : creates a vcf file from these positions