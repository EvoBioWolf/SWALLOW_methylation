### Genetic vs. Epigenetic Comparison (genome-wide)

This folder contains a suite of RMarkdown scripts designed to integrate epigenetic differentiation ($P_{ST}$) with population genetic metrics ($F_{ST}$, $\pi$, and $D_{XY}$) to investigate the drivers of methylation divergence in the Barn Swallow.

---

## Overview

### 1. CpG Island & Flanking Regions

**`1_pst_genetic_comparison_cpgi.Rmd`**

- **Focus:** Regional analysis of defined CpG Islands (CPGi) and Shores (flanking regions).
  
- **Analysis:** Compares $P_{ST}$ values against:
  
  - **$F_{ST}$**: Relative genetic differentiation.
    
  - **$\pi$ (Nucleotide Diversity)**: Within-population genetic variation.
    
- **Goal:** To determine if epigenetic divergence in regulatory regions is coupled with localized genetic structure.
  

### 2. Genome-Wide (Site-Level) Analysis

**`1_pst_genetic_dxy_comparison_genomewide.Rmd`**

- **Focus:** The comprehensive CpG site dataset (unfiltered by regional annotation).
  
- **Analysis:** Compares site-level $P_{ST}$ to genetic metrics.
  
- **Goal:** Provides a global view of the genetic-epigenetic relationship without grouping into islands/flanking regions.
  

### 3. Absolute Divergence ($D_{XY}$) Comparisons

**`1_pst_dxy_genetic_comparison_cpgi.Rmd`**

- **Focus:** Comparing $P_{ST}$ specifically to **$D_{XY}$** (Absolute Genetic Differentiation).
  
- **Rationale:** Unlike $F_{ST}$, $D_{XY}$ is not influenced by within-population diversity, making it a more robust measure of the time since divergence or ancestral polymorphism.
  
- **Analysis:** These scripts assess whether regions of high absolute genetic divergence also exhibit high methylation divergence, identifying potential regions of differentiation that are consistent across both markers (in our case not found). This is applied to the genomewide dataset in the script above.
  
### 4. Plotting

**`2_pst_rho_plots`** & **`2_pst_dxy_rho_plots`**

- **Focus:** Plotting DXY and FST correlations according to category and the number of sites represented. 
  
- **Rationale:** Because weaker correlations may be due to lower number of sites included, I want to plot correlations per category where it is clear the proportion of sites represented, in order to visualize if there is a relationship.
  
- **Analysis:** These scripts plot the correlation per category (TTS, promoter, etc.) and also include a size compnent to look at number of site represented.
---

### Prerequisites

Before running these scripts, ensure the following data objects are available:

1. **PST Tables**: Calculated via the `Pstat` package (found in the `2_Divergence/1_PST/` subfolders).
  
2. **Genetic Metrics**: $F_{ST}$, $\pi$, and $D_{XY}$ calculated from the VCF/genotype data.
  
3. **Homer Annotations**: Required for the `cpgi` scripts to link sites to specific island IDs.