3_ld_decay: Genomic LD Preparation

This directory contains the pipeline for converting genomic variant calls (VCF) into PLINK formats and calculating the baseline genetic Linkage Disequilibrium (LD) decay for the population.

### 1. Format Conversion & Pre-processing

The initial step involves transitioning from VCF files to the map/ped format required by PLINK for LD calculations.

- **`0_vcf_2_map.sh`**: The primary shell script that handles the conversion of VCF files into PLINK-compatible `.map` and `.ped` files.
  
- **`0_HR_168_plink.ped/map`**: The resulting PLINK files containing genetic information for 168 individuals.
  
- **`1_plink.keep.txt`**: A filtering file used to subset specific individuals or populations for the LD analysis.
  

### 2. LD Calculation

- **`1_plink.ld.sh`**: Executes the PLINK LD analysis. This script calculates the squared correlation coefficient ($r^2$) between pairs of SNPs within defined physical windows.
  
- **Parameters**: Typically utilizes specific window sizes and $r^2$ thresholds to capture the rate of linkage decay across the genome.
  

### 3. Result Binning

- **`2_ld_decay.py`**: A Python script used to aggregate the raw PLINK LD results.
  
- **Method**: It organizes the $r^2$ values into distance-based bins (e.g., 100bp or 1kb increments) to facilitate the plotting of a smooth decay curve.
  
- **Note**: This specific output was used for exploratory visualization and was not included in the final manuscript analysis.