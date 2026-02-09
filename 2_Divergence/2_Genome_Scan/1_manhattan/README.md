# Manhattan Plot Visualizations (`1_manhattan`)

This folder contains scripts for generating genome-wide Manhattan plots to visualize and compare epigenetic ($P_{ST}$) and genetic ($F_{ST}$) differentiation outliers.

## Overview

The primary goal of these visualizations is to identify differntiation between PST and FST. By plotting $P_{ST}$ and $F_{ST}$ across the chromosomes, we can visually assess whether epigenetic outliers co-localize with known genetic outliers (top 1% of FST distribution).

## Script Documentation

### `1_manhattan_plot_PST_sites.Rmd`

This is the core script for visualizing site-level differentiation.

- **Input**: Site-level $P_{ST}$ values (calculated via `Pstat`) and $F_{ST}$ values.
  
- **Functionality**:
  
  - Generates genome-wide Manhattan plots with chromosomes ordered numerically.
    
  - Highlights **Outliers**: Typically defined as the top 1% of the distribution.
    
  - **Comparison**: Allows for a vertical comparison between genetic ($F_{ST}$) and epigenetic ($P_{ST}$) peaks to identify regions of synchronized divergence.