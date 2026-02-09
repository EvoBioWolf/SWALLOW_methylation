# 2_PST_nearest_FST

This folder contains information on determining if high PST windows/sutes lie near FST peaks. I want to determine if epigenetic outliers are closer to genetic outliers more often than expected by chance.

## Workflow

### 1. Analysis & Simulation

- **`1_PST_near_FST_simulation_*.r`**
  
  The core scripts of this analysis. Here, either for sites, cpgi, or windows I calculate the genomic distances between $P_{ST}$ and $F_{ST}$ peaks of the empirical dataset. I then perform permutation tests using randomized $P_{ST}$ datasets to generate a null distribution, allowing for the calculation of significance for the observed distances. I always compare the top 1% of outliers in FST and PST.
  
- **`./2_simulation_results/`**
  
  Stores the output of the simulation script.
  
  - Files prefixed with `1_PST...` contain the results of the randomized simulations.
    

### 2. Visualization

- **`2_plot_PST_nearest_FST_random_PST_gw.r`**
  
  Generates visualizations for the observed vs. simulated distributions for the Figures. Here, you can change the script from the simulatons folder to look at output for different runs.