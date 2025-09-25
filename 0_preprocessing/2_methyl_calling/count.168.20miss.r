###09.08.22
##Sarah Mueller
##script to make meth count files needed for further analysis

##directory for libraries
.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))

#Libraries
library(MASS)
library(vegan)
library(viridis)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(broom)
library(readxl)
library(cluster)
library(ape)
library(ggrepel)
library(ggpubr)
library(lme4)
library(data.table)
library(methylKit)
library(emmeans)
library(car)


setwd("/dss/dsshome1/lxc0B/ra52qed/output/1.RRBS/4.0.methylunite/")
#sample names
sample_names <- read.table("/dss/dsshome1/lxc0B/ra52qed/scripts/0.sample.lists/sample.names.168.txt")
#assign meth, total and unmethylated to sample names
sample_names <- cbind(sample_names, sample_names, sample_names)
colnames(sample_names) <- c("ME", "TO", "UN")
sample_names$ME <- sub("^", "ME_", sample_names$ME)
sample_names$UN <- sub("^", "UN_", sample_names$UN)
sample_names$TO <- sub("^", "TO_", sample_names$TO)

#sample size
N =168

##I want to have a dataframe with the sample names and methylation data, as a proportion and counts of methylation
##Durign this I will also remove low variablility sites, so the outcome of this script is what I will run further data analysis on
################################################
##load methylation data
load("meth.168.20miss.RData")
meth_basic_stat <- methylKit::getData(meth.168.20miss)
##get proportion data by taking number of methylated/total coverage
proportion_basic_stat <- (meth_basic_stat[ , c(unique(c(seq(6, ncol(meth_basic_stat), 3))))]/meth_basic_stat[ , c(unique(c(seq(5, ncol(meth_basic_stat), 3))))])
#change names to sample names for easy identification
colnames(proportion_basic_stat) <- sample_names[,1]
meth_basic_stat<- meth_basic_stat %>%
  tidyr::unite(site.id, chr:start)
rownames(proportion_basic_stat) <- meth_basic_stat[,1]

#############################
##filter out sites with low variability (<.10 and >0.9, this is sort of a standard from papers, but could be changed)
#############################
proportion_basic_stat<- proportion_basic_stat %>%
  mutate(proportion.ave = rowMeans(proportion_basic_stat, na.rm = TRUE))
proportion_basic_stat<- proportion_basic_stat %>%
  mutate(site.id = meth_basic_stat[,1])
proportion_basic_stat<- proportion_basic_stat %>%                       
  rowwise %>%
  filter(proportion.ave > 0.01) %>%
  filter(proportion.ave < 0.99)

#It is generally better to use the raw count data, along with total read count or unmethylated reads
#So what I want to do is take the sites that have been filtered by proportion and use these to extract the count data at these sites

##methylated
meth_count <- left_join(proportion_basic_stat, meth_basic_stat, by="site.id")
info_col <- meth_count[,(N+2)]
meth_count <- (meth_count[ , c(unique(c(seq((N+6), ncol(meth_count), 3))))])
colnames(meth_count) <- sample_names[,1]

##we have to extract the total counts from the datatset.
meth_total <- left_join(proportion_basic_stat, meth_basic_stat, by="site.id")
meth_total <- (meth_total[ , c(unique(c(seq((N+5), ncol(meth_total), 3))))])
colnames(meth_total) <- sample_names[,2]

##can also add unmethylated if needed
meth_count_un <- left_join(proportion_basic_stat, meth_basic_stat, by="site.id")
meth_count_un <- (meth_count_un[ , c(unique(c(seq((N+7), ncol(meth_count_un), 3))))])
colnames(meth_count_un) <- sample_names[,3]


##combine these into one dataframe
count.168.20miss <- cbind(info_col,meth_count, meth_total, meth_count_un)
##combine proportion data
prop.168.20miss <- cbind(info_col, proportion_basic_stat[1:168])

#remove sites associated with sex (see sex.evaluation folder)
sex.dmc <- read.table(file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_RRBS/3.2.sex.evaluation/167.dmc.sex.mean.bed")
#merge chr and position into one site.id so this can be used in an antijoin with the proportion and count tables
sex.dmc$site.id <- paste(sex.dmc$V1, sex.dmc$V2, sep = ".")
#use anti-join to remove sex-related sites from the dataset
prop.168.20miss.nosex <- anti_join(prop.168.20miss, sex.dmc, by="site.id")

#perform again on the count table
count.168.20miss.nosex <- anti_join(count.168.20miss, sex.dmc, by="site.id")

##save
save(count.168.20miss, file = "count.168.20miss.RData")
save(count.168.20miss.nosex, file = "count.168.20miss.nosex.RData")
save(prop.168.20miss, file = "prop.168.20miss.RData")
save(prop.168.20miss.nosex, file ="prop.168.20miss.nosex.RData")

