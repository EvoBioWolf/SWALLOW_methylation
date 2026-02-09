##Using DSS to identify age associated sites and sex associated sites
##Sarah A. Mueller
#16.09.22 started

.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))

library(DSS)
require(bsseq)
library(tidyr)
library(dplyr)
library(maditr)

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/3_sex_evaluation")

#first load in counts from the Dataset 3, because this should include all sites in Dataset1 &2 so we only have to do the analysis once to get sites to clean from all three datasets
#sample names
sample_names <- read.table("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.168.txt")

N=168

#read in methylation data (average methylation at CpG loci per individual)
prop.168.merged.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/prop.168.merged.20miss.txt", header=TRUE)
meth.168.merged.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/meth.168.merged.20miss.rmlowv.rmout.txt", header =TRUE)

##methylated
meth_count <- (meth.168.merged.20miss[ , c(unique(c(seq((5), ncol(meth.168.merged.20miss), 3))))])
colnames(meth_count) <- sample_names[,1]
meth_count <- meth_count %>% 
  dplyr::mutate(site.id = paste(meth.168.merged.20miss$site.id)) 
meth_count_long <- meth_count %>% pivot_longer(cols=!site.id, values_to = "X", names_to = "ID")

#total
meth_total <- (meth.168.merged.20miss[ , c(unique(c(seq((4), ncol(meth.168.merged.20miss), 3))))])
colnames(meth_total) <- sample_names[,1]
meth_total <- meth_total %>% 
  dplyr::mutate(site.id = paste(meth.168.merged.20miss$site.id)) 
meth_total_long <- meth_total %>% pivot_longer(cols=!site.id, values_to = "N", names_to = "ID")
#merge
meth <- full_join(meth_count_long, meth_total_long, by=c("site.id", "ID"))
###

#make chr and pos two columns
dss_input <- meth %>% 
              tidyr::separate(site.id, c('chr', 'pos'))
#make sure start and end are treated as integers
dss_input$pos <- as.integer(dss_input$pos)
#keep only columns we will work with
dss_input <- dss_input %>%
  dplyr::select("chr","pos", "N", "X", "ID")

#remove NA
dss_input_clean <- dss_input %>%
  filter(!is.na(X) & !is.na(N))

#########################
###############PROBLEM WITH MISSING DATA

#make allnames 
allnames <- NULL

#loop through them..
for (samp in unique(dss_input_clean$ID)) {
  
  var <- dss_input_clean[grepl(samp,dss_input_clean$ID),][,1:5]
  
  cat('Saved to variable: ',samp,'\n')
  #assign it to a new variable
  assign(paste0(samp),var)
  allnames <- unique(c(allnames,samp))
  
}

#get names
vars <- mget(allnames)

#create object
BSobj = makeBSseqData(vars,c(allnames))
BSobj

save(BSobj, file = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/3_sex_evaluation/2_168_sex_eval_BSobj.RData")
