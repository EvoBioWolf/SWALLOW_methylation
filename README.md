# Covariation of genetic and epigenetic divergence in barn swallows

## Authors

Sarah A. Mueller, Uthara Srinivasan, Zachary M. Laubach, Drew R. Schield, Jan A. C. von Rönn, Rebecca J. Safran, Jochen B. W. Wolf

## Project Description

This repository contains the analysis pipeline and custom scripts for the study titled **"Covariation of genetic and epigenetic divergence in barn swallows"**.

The project investigates the population genomic structure and epigenetic patterns of Barn swallow (Hirundo rustica) using Reduced Representation Bisulfite Sequencing (RRBS) and Whole Genome Sequencing (WGS) data. Our workflow covers the entire process from raw sequence pre-processing and methyl-calling to population structuring analyses including PCA, Redundancy Analysis (RDA), $P_{ST}$ (Phenotypic Fixation Index) calculations, linkage evaluation and hybrid zone analysis.

---

## Repository Structure

The repository is organized following the sequential steps of the analysis:

- **`0_Pre-processing/`**: Scripts for adapter trimming, quality control, and methyl-calling.
  
- **`0_genome_files/`**: Reference genome metadata and filtering scripts (e.g., C-T polymorphism filters).
  
- **`0_metadata/`**: Sample manifests, mapping coordinates, and environmental data.
  
- **`1_Population_Structure/`**: Analysis of genetic and epigenetic clustering (PCA, Admixture).
  
- **`2_Divergence/`**: Calculations of $F_{ST}$ and $P_{ST}$ to identify genomic regions under selection.
  
- **`3_Linkage/`**: Linkage Disequilibrium (LD) decay and haplotype analysis.
  
- **`4_RGT_Hybrids/`**: Specific scripts for hybrid zone analysis and introgression.
  
- **`5_WGS/`**: Comparative analysis using Whole Genome Sequencing data.
  

---

## Data Availability

Due to the large size of genomic data files (e.g., `.ped`, `.vcf`, and `.bam` files), the raw data is not hosted directly in this GitHub repository.

- **Raw Data:** Available via NCBI BioProject under accession number PRJNA797941.
  

---

## Contact

For questions regarding the code or data, please open an **Issue** in this repository or contact **Sarah Mueller** at `s.mueller@bio.lmu.de`.