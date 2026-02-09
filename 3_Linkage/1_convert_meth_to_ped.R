### Libraries
.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))

library(tidyverse)
library(ggplot2, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")

####Load Datasets######

prop.168.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/prop.168.merged.20miss.rmlowv.rmout.txt", header=T)

# Convert numeric columns based on the given conditions
HR.168.ped <- prop.168.20miss %>%
  mutate(across(where(is.numeric), 
                ~ case_when(
                  is.na(.) ~ "00", 
                  . < 0.3 ~ "AA", 
                  . >= 0.3 & . <= 0.7 ~ "AG", 
                  . > 0.7 ~ "GG"
                ), 
                .names = "{col}_converted"))  # Creates new columns with the '_converted' suffix

HR.168.ped <- HR.168.ped %>%
        select("site.id", ends_with("_converted")) %>%
        select("site.id", starts_with("ME_"))

# Remove prefix 'ME_' and suffix '_converted' from column names
HR.168.ped <- HR.168.ped %>%
  rename_with(~ gsub("^ME_|_converted$", "", .))

HR.168.t.ped <- as.data.frame(t(HR.168.ped))

HR.168.t.ped <- HR.168.t.ped %>%
  mutate(
    `Individual ID` = rownames(HR.168.t.ped),  # Set Individual ID as the row names
    `Paternal ID` = 0,                         # Set Paternal ID to 0
    `Maternal ID` = 0,                         # Set Maternal ID to 0
    Sex = "unknown",                           # Set Sex to "unknown"
    Phenotype = 0                              # Set Phenotype to 0
  ) %>%
  select(`Individual ID`, `Paternal ID`, `Maternal ID`, Sex, Phenotype, everything())  # Move mutated columns to the front


#Remove first row which contains site.id (this will be in .map file below)
HR.168.t.ped <- HR.168.t.ped[-1,]

HR.168.map <- HR.168.ped %>%
  select(site.id) %>%
  mutate(
    chromosome = sapply(strsplit(as.character(site.id), "_"), `[`, 1),
    snp = sapply(strsplit(as.character(site.id), "_"), `[`, 2),
    Genetic_distance = 0,  # Set genetic distance to 0
    Base_pair_position = snp  # Base-pair position is same as SNP
  ) %>%
  select(chromosome, snp, Genetic_distance, Base_pair_position)


write.table(HR.168.t.ped, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.ped", row.names=F, col.names=F, quote=F)
write.table(HR.168.map, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/0_HR.85.20miss.rmlowv.rmout.map", row.names=F, col.names=F, quote=F)
