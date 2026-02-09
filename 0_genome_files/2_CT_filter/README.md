## SNP .bed file for filtering

I wanted to pull all SNP loci that have a C or G in the reference genome to remove any SNPs that would interfere with the CpG site.

CpG - if either the C or the G is replaced by a SNP, this will cause those with the SNP variant to be unable to be methylated. This can cause spurious associations of 'methylation' to the genomic background, when it is really the signal of a SNP loci.

To be conservative, I pulled variants from both the vcf file created from all 168 individuals, as well as the vcf file from Schield et al. 2024, to catch variants that may not have been sequenced in the European populations and would therefore not be included in my vcf file, but may still be present in the other populations.

Run 1_extract_all_snps_HR.168.sh to get a bed file that is used for downstream filtering.