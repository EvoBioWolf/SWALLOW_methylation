1_chr11: Inversion LD Analysis

This sub-folder focuses on the structural variant on Chromosome 11. Because inversions act as recombination suppressors, we expect to see distinct patterns of Linkage Disequilibrium (LD) and associated methylation "haplotypes" in this region.

### 1. Data Subsetting & PLINK LD

The analysis begins by isolating the Chromosome 11 inversion region and calculating LD specifically for this structural variant.

- **`0_thin_SNPs.sh`**: thins the SNP density to reduce computational load and minimize redundant signal for the linkage analysis.
  
- **`1_plink_ld`**: Runs the PLINK LD command on the inversion region coordinates from the other folder.
  
- **`HR.85.maf2...` files**: The resulting PLINK binary files (`.bed`, `.bim`, `.fam`) and LD output (`.ld.gz`) for the 85-individual cohort, filtered for MAF > 0.2.
  

### 2. Result Partitioning

- **`1_separate_results.sh`**: Separates the bulk PLINK LD output into three categories:
  
  - **`HR_85_snp...`**: SNP-SNP linkage.
    
  - **`HR_85_meth...`**: CpG-CpG (methylation) linkage.
    
  - **`HR_85_snp.meth...`**: Cross-variant (SNP-CpG) linkage.
    

### 3. Visualization & Heatmaps

To visualize the block structure of the inversion, we generate triangular linkage heatmaps.

- **`2_all_linkage_heatmaps.r`**: The R script used to generate LD heatmaps.
  
- **`2_all_chr11_inv_heatmap_ordered.png`**: The primary heatmap output showing the structured blocks of LD within the inversion.
  
- **`3_high_LD_haplotypes[_label].png`**: Visual representations of the specific high-LD blocks (haplotypes) identified within the inversion region.
  

### 4. Sliding Window Analysis

To observe how $r^2$ changes across the chromosome, data is aggregated into 10kb windows.

- **`6_create_10kb_r2_windows.r`**: Processes the LD data into spatial windows.
  
- **`6_..._10kb_windows.txt`**: The resulting windowed data for SNPs, Methylation, and SNP-Meth pairs, allowing for the plotting of $r^2$ values relative to the physical position on Chromosome 11.