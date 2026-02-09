3_Linkage

In this folder, I focus on quantifying the linkage between genetic variants (SNPs) and epigenetic variants (CpG methylation), as well as between CpG sites and between SNPs. By converting methylation levels into a pseudo-genotype format, we can calculate Linkage Disequilibrium (LD) across the genome.

## 1. Data Conversion & Formatting

- **`1_convert_meth_to_ped.R`**: This script transforms methylation proportion data into a categorical "genotype" format suitable for PLINK. It generates the `.map` and `.ped` files for the methylation sites. The primary goal is to transform methylation data into PLINK-compatible formats (`.map` and `.ped`) to merge them with existing Whole Genome Sequencing (WGS) data.
  
- **`0_HR.85.20miss...` files**: These represent the base genetic and combined SNP-methylation datasets for the 85-individual cohort.
  
- To merge the converted methylation files with the existing genetic data (referenced from the WGS project, look at scripts within: `/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/3_ld_decay`).
  
- The result of merging scripts are the files beginning with `1_HR.85.20miss` which are the combined SNP and CpG sites.
  
- `run_r.sh`: A wrapper script to execute the R-based portions of the pipeline.
  

### 2. Linkage Disequilibrium (LD) Calculation

- **`1_plink.ld.sh`**: Runs the LD calculations using PLINK.
  
  - **Filtering**: We apply a Minimum Allele Frequency (**MAF**) filter of **> 0.2**. This high threshold ensures a sufficient number of individuals are present in each "allele" group (methylation state) to make robust statistical inferences regarding linkage.