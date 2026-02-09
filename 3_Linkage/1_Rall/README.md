## 1_Rall: Partitioned LD Analysis

This folder contains the secondary analysis of Linkage Disequilibrium (LD) results. It takes the global output from PLINK and partitions it to compare the decay patterns of genetic and epigenetic variants, specifically across different genomic regions (e.g., islands, promoters, exons).

### 1. Data Partitioning

The raw output from the PLINK LD analysis (the `HR_Rall` compressed files) is separated into three distinct biological categories to compare how linkage differs between variant types.

- **`1_separate_results.sh`**: The script used to parse the raw PLINK output into the categorized files below.
  
- **`2_Rall_snp_maf2...txt`**: LD results for **SNP-SNP** pairs.
  
- **`2_Rall_meth_maf2...txt`**: LD results for **CpG-CpG** (methylation) pairs.
  
- **`2_Rall_snp_meth_maf2...txt`**: LD results for **SNP-CpG** cross-linkage.
  

### 2. Plotting and Exploration

Once partitioned, the data is analyzed to visualize LD decay—the rate at which $r^2$ decreases as physical distance between markers increases.

- **`2_ld_decay_plotting.r`**: Generates the primary LD decay curves for the three categories.
  
- **`2_Rall_LD_exploration.[Rmd/html]`**: An interactive analysis notebook used to explore the specific distributions of $r^2$ values and validate the results.
  

### 3. Regional LD Analysis

We extend the analysis beyond global patterns to see how LD is influenced by genomic architecture and recombination rates.

- **`3_all_LD_regions.Rmd`**: Analyzes LD patterns across specific genomic contexts:
  
  - **Methylation context**: CpG Islands, shores/flanks, and open sea sites.
    
  - **Gene architecture**: Promoters, Transcription Termination Sites (TTS), Exons, Introns, and Intergenic regions.
    
- **`rustica.recomb.rename.chr.txt` / `rustica.rmap.100kb.txt`**: Recombination maps used to calibrate LD expectations against known recombination hotspots and coldspots in the *rustica* population.