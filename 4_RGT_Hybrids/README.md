# 4_RGT_Hybrids: Hybrid Methylation Analysis

This directory contains the workflow for investigating epigenetic patterns in hybrid individuals compared to their parental populations. It integrates genomic hybrid indices and local ancestry to model the drivers of methylation.

### 1_PCA: Population Clustering

- **Purpose**: Broad visualization of the 144 individuals, including all three parental populations and their respective hybrids using prcomp in R.
  

### 2_Parental_high_PST: Outlier Visualization

This folder focuses on the top 1% of $P_{ST}$ outliers identified across different genomic features (Islands, Flanks, and Sites).

- **`1_RG_PST_heatmaps_nohybrids.[Rmd/html]`**: Heatmaps for the RG hybrid zone parental populations.
  
- **`1_RT_PST_heatmaps_nohybrids.[Rmd/html]`**: Heatmaps for the RT hybrid zone parental populations.
  
- **Visuals**: `Figure3_RG_Heatmap.png` and `Figure3_RT_Heatmap.png` show the average methylation differences that define the parental divergence before assessing hybrid states.
  

### 3_methylation_models: Architecture & Modeling

This is the core analytical folder, subdivided by hybrid zone (**RG_hyb** and **RT_hyb**).

#### A. Differential & Inheritance Analysis (`RG_diff` / `RT_diff`)

- **Hybrid Classification**: Uses classes (e.g., F1, F2, Backcross) identified via *NewHybrids* (see WGS project directory).
  
- **Statistical Testing**: Runs **ANOVA** with post-hoc **Tukey** tests to compare parentals and F1 hybrids.
  
- **Genomic Architecture**: Determines inheritance patterns:
  
  - **Dominance**: Hybrid methylation matches one parent.
    
  - **Co-dominant**: Hybrid methylation is intermediate.
    
  - **Transgressive**: Hybrid methylation is outside the range of both parents.
    
- **Files**: `0_RG_parental_hybrid_tukey_results.txt` and associated boxplots.
  

#### B. Prediction Models (`RG_models` / `RT_models`)

- **Purpose**: For significantly different sites, this analysis determines the best predictor of methylation values in hybrids.
  
- **Predictors tested**:
  
  1. **Nearest SNP**: Direct genetic cis-influence.
    
  2. **Local Genomic Context**: Local ancestry (ancestry of the specific genomic block).
    
  3. **Global Hybrid Index**: Overall genomic composition (derived from *introgress*).
    
- **Key Files**:
  
  - `1_RG_models.Rmd`: The main modeling script.
    
  - `2_Figure3e_local_anc_plot.png`: Visualization of local vs. global ancestry effects.
    

---

## 🧬 Data Dependencies

This workflow relies on genetic data generated in the WGS pipeline:

- **Hybrid Indices**: `/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids`
  
- **Local Ancestry/Introgress**: `/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/1_introgress`