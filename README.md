# Genetic and Epigenetic variation in barn swallows

Mueller et al. XXXX

This Github repository contains scripts for the article titled "Covariation of genetic and epigenetic divergence in barn swallows". A small overview is below:

#### 0_genome_files

This folder contains information on how to create files related to the genome of the barn swallow inclusing cut site filtering, CpG islands identification. Most of these files are then used for pre-processing and filtering of methylation calls.

#### 0_preprocessing

Contains information on bioinformatic processing of RRBS data from fastq files to the final datasets used in the manuscript for analysis. For further details see .README inside this folder.

#### 1_divergence

Contains the script to calculate site/region PST between populations using Pstat

#### 2_linkage

Contains scripts for converting methylation into an 'epi'-genotype (0,1,2) so that LD analysis between SNPs and CpG sites can be carried out

#### 3_hybrids

Contains scripts for analysing individuals from hyrbid zones to determine significantly different methylation regions/sites between parental populations and quantify the impact of cis- and trans- contribution.
