NewHybrids Analysis

I used the software **NewHybrids** to identify hybrid individuals and their specific hybrid classes (e.g., F1, F2, or Backcrosses) within the swallow population. By analyzing genetic markers (Ancestry Informative Loci), we can distinguish between pure parental species and different generations of hybridization. This is done twice: 1 run for each of the hybrid zones (rustica/gutturalis and rustica/tytleri)

#### 1. Empirical Data Processing (`1_newhybrids_script.R`)

This script prepares the raw genomic data (VCF files) for analysis.

- **Format Conversion:** It converts genetic genotypes into the specific numerical format required by NewHybrids.
  
- **Execution:** It runs hybriddetector and NewHybrids algorithm on the simulated and empirical samples collected from the field.
  

#### 2. Simulation and Validation (`2_newhybrids_results_sim.r`)

- **Performance Testing:** We test how well the software identifies the known individuals from simulations. This allows us to calculate **Efficiency** (how many hybrids we caught) and **Accuracy** (how many we identified correctly).
  

- **Threshold Selection:** We test multiple probability cut-offs (from 0.8 to 0.995) to find the "sweet spot" that gives us the most reliable results for our specific dataset.
  

### Key Outputs

- **Performance Plots:** Visualizations showing how well the software performs across different hybrid classes.
  
- **Final Assignments:** A table (`2_output_classes.txt`) assigning each swallow to a genetic category (Pure Parent, F1, F2, etc.) based on the optimized probability thresholds.