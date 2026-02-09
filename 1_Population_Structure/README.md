# Population Structuring Workflow

This directory contains information for the Population Structuring analyses. Here, I use PCA, RDA and Mantel tests to look at overall population structuring in methylation data. For equivalent analyses in genetic data see `/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS` .

## 📂 Directory Structure Overview

| **Folder** | **Description** |
| --- | --- |
| `1_PCA` | Principal Component Analyses (Unsupervised clustering) |
| `2_RDA` | Redundancy Analyses (Environment-constrained ordination) |
| `3_Mantel` | Mantel tests and partial Mantel tests |

---

## 1_PCA: Principal Component Analysis

This folder contains scripts and outputs for unsupervised clustering based on methylation proportions using `prcomp`.

- **1_PCA_85_merged_20miss_nolowv_rmout.[Rmd/html]**: The primary analysis for all 85 individuals from pure subspecies, excluding sites with missing data. This can be adapted for running in specific regions as well (only islands, flanks, etc.)
  
- **1_PCA_chr11_LD_85_merged_20miss_nolowv_rmout.[Rmd/html]**: PCA focused specifically on **Chromosome 11**, targeting a known inversion to observe localized methylation structuring.
  
- **1_PCA_rustica_merged_20miss_nolowv_rmout.[Rmd/html]**: Analysis restricted to the **rustica** population only, used to determine sub-structuring for supplementary figures.
  

## 2_RDA: Redundancy Analysis

This folder focuses on constrained ordination to explain methylation variance using environmental predictors. Overall, we use proportion methylation as response and 8 geographic, genetic and environmental predictors.

- **2_RDA.85.20miss.RRBS_cpgi.[Rmd/html]**: RDA analysis in methylation within cpgialsnds (done to see sensistivity of results between regions)
  
- **2_RDA.85.20miss.RRBS.Rmd/html**: Core analysis done to look at methylation structuring in 85 individuals, across all sites that do not contain missing data. This was performed because PC1 in the unsupervised PCA was dominated by noise.
  

## 3_Mantel: Mantel Tests & Data Prep

This folder contains the preparation of genetic objects and the subsequent correlation tests between distance matrices.

### Data Preparation (Genlight/Genepop)

- **0_genlight_obj.r**: The initial script to create a `genlight` object from raw genetic data.
  
- **0_genlight.85.10miss.RData**: The resulting RData object for individuals (10% missingness threshold).
  
- **0_genpop_85_10miss.RData**: The genetic data converted to `genepop` format for population-level analysis.
  

### Statistical Testing

- **1_mantel_meth_distance.Rmd**: Script for calculating Manhattan and Euclidean distances for methylation data.
  
- **d.gen.ind.euc.10miss.RData**: Euclidean genetic distance matrix (Individual level).
  
- **d.gen.pop.manual.168.5miss.RData**: Genetic distance matrix (Population level).
  
- **run.r.sh**: Shell script to automate the R processing.