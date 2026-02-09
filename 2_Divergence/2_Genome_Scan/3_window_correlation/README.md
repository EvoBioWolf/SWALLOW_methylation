### Correlation and PST Analysis Documentation

Here, I summarize the scripts used for testing genomic window parameters and evaluating the relationship between  $P_{ST}$  and $F_{ST}$.

---

### 0. Parameter Testing

### `0_correlation_window_parameter_testing.Rmd`

This script is used for testing different window sizes (e.g., 1kb, 5kb, 10kb) to assess their impact on correlation estimates between $P_{ST}$ and **$F_{ST}$**.

- **Key Finding:** Generally, lower window sizes lead to a slight drop in correlation, though the overall trend remains consistent.
  
- **Caveat:** The observed lower correlation at smaller scales may be a result of not re-calculating $F_{ST}$ specifically for those smaller window sizes.
  

---

### 1. Methylation Distance Correlation

### `1_correlation_10kb_meth.dist.Rmd`

Calculates correlations between $F_{ST}$ and $P_{ST}$ for **10kb windows** across various genomic contexts.

#### Genomic Contexts Analyzed:

- **CpG Context:** Islands, shores, and openseas.
  
- **Functional Annotations:** Sub-categories include intron, exon, promoter, TTS, and intergenic regions (for both CpG and non-CpG regions).
  

#### Prerequisites:

1. **Regional Estimates:** Must first be calculated using:
  
  `/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/1_CpG_islands/1_regional_windows.Rmd`
  
2. **PST Calculations:** PST for each window must be calculated beforehand. Scripts are located in:
  
  `/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_window/`
  

---

## 2. Distributional Comparisons

### `2_PST_FST_distribution.Rmd`

This script is not included in the final paper but serves as an exploratory analysis of how $P_{ST}$ values shift in relation to $F_{ST}$ levels.

- **Analysis:** Compares the distribution of $P_{ST}$ when $F_{ST}$ is high versus background $F_{ST}$ levels.
  
- **Goal:** To identify potential shifts in the $P_{ST}$ distribution triggered by shifts in $F_{ST}$.