.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))

library(vegan)
library(tidyverse)
library(readxl)
library(cluster)
library(vcfR)
library(adegenet)
library(SNPRelate)
library(gdata)
library(dartR)

#sample names
sample_names <- read.table("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.85.sub.pop.txt", header=T)
sample_names_v <- sample_names$V1
#get explanatory variables set
info <- read_xlsx("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/NCBI_info_updated_20220727_climate_data.xlsx")
#remove all individuals not used for this analysis
info=info %>%
  dplyr::slice(match(sample_names$V1, info$sample_name))

##genetic matrix
vcf <- read.vcfR("/dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/HR.85.biallelic.maf05.dp_3_60.vcf.gz")
#convert to gen_light object
gen_light <- vcfR2genlight(vcf)

#add pop information
pop(gen_light) <- as.factor(c(sample_names$sub))
popNames(gen_light)
ploidy(gen_light) <- 2

#save for future use
save(gen_light, file= "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/1_Population_Structure/3_Mantel/0_genlight.85.10miss.RData")
