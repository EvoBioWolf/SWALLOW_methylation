# 1_Population_Structuring (Genetic Data)

This directory contains the genetic counterpart to the methylation structuring analysis. It utilizes Whole Genome Sequencing (WGS) data to characterize population differentiation ($F_{ST}$), absolute divergence ($D_{XY}$), and nucleotide diversity ($\pi$).

### 1. Data Cleaning & Neutral Structure

- **`0_LD_prunning/`**: Contains scripts to prune SNPs based on Linkage Disequilibrium. This ensures that the PCA reflects broad population structure rather than localized linkage blocks (like inversions).
  
- **`1_PCA/`**: Principal Component Analysis of the genetic data to visualize the major axes of variation among the 168 individuals.
  
- **`1_168_pop_map.txt`**: The master population map linking individual IDs to their respective subspecies and geographic locations.
  
- **`3_FST/`**: Relative differentiation calculated between population pairs.
  

### 2. Diversification Statistics (pixy)

We use `pixy` to calculate several key population genetic statistics simultaneously. This is preferred over standard methods as it correctly handles invariant sites.

- **`4_pixy/`**: The raw output folder for the pixy analysis.
  
  - **`2_PI/`**: Nucleotide diversity ($\pi$) results, quantifying the genetic health and effective population size of each group.
    
  - **`2_FST/`** Relative differentiation calculated between population pairs.
    
  - **`2_DXY/`**: Absolute genetic divergence, used to identify regions of the genome that diverged before population splitting.
    
  - **`2_pixy_plot.Rmd` / `2_pixy_plot_no_pi.Rmd`**: RMarkdown scripts for generating Manhattan plots of $F_{ST}$, $D_{XY}$, and $\pi$ across the genome.
    
  - **`1_pistats.sh`**: Shell script to automate the extraction of summary statistics from the diversity outputs.
    
  - **`3_highdxy_high_fst.r`**: A specialized script to identify "genomic islands of divergence" where both relative and absolute differentiation are elevated.