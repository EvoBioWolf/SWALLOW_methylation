library(dplyr)
library(tidyr)
library(Pstat)
library(readxl)
library(zoo)
library(vegan)

#______Sample Information______
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_cpgi_pst/cpgi")
#sample names
sample_names <- read.table("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.85.txt")
#get explanatory variables set
info <- read_xlsx("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/NCBI_info_updated_20220727.xlsx")

#_________All Sites Dataset________
#load all-site methylation data
prop.20miss <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/1_CpG_islands/168.20miss/homer.cpgi.prop.168.merged.20miss.rmlowv.rmout.txt", header=T)

#do not keep hybrids
prop.20miss <- prop.20miss %>%
  # 1. Move site.id to the first position
  dplyr::relocate(region.id) %>%
  # 2. Remove columns containing any of your target strings
  dplyr::select(-contains("_RG_"), 
                -contains("_RT_"), 
                -contains("proportion"),
                -contains("cpgs"))

#dataframe must be a row for each individual, with the first column containing population information and the rest containing trait variables (in this case epigenetic data)
p.df = t(prop.20miss[,2:86])

#dataframe
p.df <- as.data.frame(p.df)
names(p.df) = prop.20miss[, 1]

##add population column to sample_names
sample_names <- sample_names %>%
  mutate(pop = case_when(
    grepl("RU_13_R_", V1) ~ "rustica.r",
    grepl("CH_", V1) ~ "rustica.eu",
    grepl("DE_", V1) ~ "rustica.eu",
    grepl("FI", V1) ~ "rustica.eu",
    grepl("_G_", V1) ~ "gutturalis",
    grepl("_E_", V1) ~ "erythrogaster",
    grepl("_S_", V1) ~ "savignii",
    grepl("_TV_", V1) ~ "transitiva",
    grepl("_T_", V1) ~ "tytleri"))

p.df <- cbind(sample_names$pop, p.df)
colnames(p.df)[1] <- "pop"

##Create a list of comparisons for PST calculation
pop.list <- as.data.frame(unique(as.vector(p.df$pop)))
colnames(pop.list) <- c("pop1")
pop.combo <- tidyr::crossing(var1=pop.list$pop1, var2=pop.list$pop1)
pop.combo <- pop.combo[c(2:7, 10:14, 18:21,26:28, 34:35,42),]

Dataset = p.df

for(i in 1:nrow(pop.combo)) {
  row <- pop.combo[i,]
  
  # Define populations to analyse
  Pop1<- row$var1
  Pop2<- row$var2
  Pops <- paste0(row$var1, ".vs.", row$var2, ".pst")
  Output_filename<- paste0("pst.allsites.", row$var1, ".vs.", row$var2, ".txt")
  
  #extract just the popultion names
  pop<-Dataset[,1]
  
  # create a vector with n elements
  row_num <- c(2:(ncol(Dataset)-1))
  
  # specify the chunk length
  chunklength<- 500
  
  # Step 1: Split the tables
  # Split up a vector with no. of elemets = number of sites (no. of columns) in barn swallow genome
  
  split_p.df <- print(split(row_num,ceiling(seq_along(row_num) / chunklength)))
  
  # Step 2: Perform the analysis on each split table
  output_tables <- list()  # Initialize an empty list to store the output tables
  
  for (i in 1:length(split_p.df)) {
    #Atribute the split up numbers to columns/ marekrs and cbind with the population names
    current_table = cbind(pop,p.df[,as.vector(split_p.df[[i]])])
    
    # Apply PST analysis code to the current table
    output_pst<-Pst(current_table, Pw=c(Pop1, Pop2), csh = 0.25)
    
    # Store the result in the output tables list
    output_tables[[i]] <- output_pst
  }
  
  # Step 3: Row bind the output tables
  combined_table <- do.call(rbind, output_tables)
  colnames(combined_table) <- c("site.id", Pops)
  
  # Step 4: Save combined
  write.table(combined_table, file= Output_filename, quote = F, sep = " ", row.names = F)
  
}
