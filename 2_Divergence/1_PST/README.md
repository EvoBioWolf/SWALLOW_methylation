# PST Analysis

This directory contains the pipeline for calculating and exploring epigenetic differentiation ($P_{ST}$) across various genomic scales in the Barn Swallow (*Hirundo rustica*).

## Directory Structure

### 1. Sub-folders: Chunked Processing

These folders contain scripts to run the `Pstat` R package. To handle large datasets, analysis is performed in parallelized "chunks."

- **`cpgi/`**: $P_{ST}$ calculations specifically for CpG Island datasets.
  
- **`sites/`**: Individual CpG-site resolution $P_{ST}$ analysis.
  
- **`windows/`**: $P_{ST}$ calculated across 10kb tiling windows.
  

### 2. Scripts: Results Exploration & Visualization

These scripts aggregate the chunked results to perform population-level comparisons and generate publication-ready figures.

**`0_PST_results_gw_comparison.Rmd`**

- **Purpose**: Genome-wide population differentiation comparison using heatmaps.
  
- **Key Functions**:
  
  - Calculates mean P_{ST} across all loci for distinct genomic contexts: **CpG Islands, Shores, and Open Seas**.
    
  - Generates symmetric **Heatmaps** for pairwise population comparisons.
    
  - **Outputs**: Publication figures (e.g., `Figure1d_PST_heatmap.png`) and Supplemental figures for specific contexts.
    

**`1_PST_results_exploration.Rmd`**

- **Purpose**: Exploratory analysis of differentiation patterns across hybrid zones.
  
- **Key Functions**:
  
  - Identifies the **Top 1%** of diverged sites for each pairwise population comparison.
    
  - Generates **Venn Diagrams** to identify overlapping highly diverged sites across different comparisons.
    
  - Calculates "overlap counts" to find candidate islands/sites appearing in >5 comparisons.
    
  - Integrates with **HOMER** annotations to determine genomic categories (Exons, Promoters, etc.) and distance to TSS.
    
  - Produces **Volcano Plots** (Methylation Difference vs. $P_{ST}$) for specific hybrid zones ($R$ vs $G$, $R$ vs $T$).