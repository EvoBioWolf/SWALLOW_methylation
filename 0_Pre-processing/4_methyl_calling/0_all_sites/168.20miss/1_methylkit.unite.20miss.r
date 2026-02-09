library(methylKit)
library(readxl)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_output_files_test")
list.files()

#CHANGE FOR EACH NEW METHYLOBJECT
#name of how you want the data to be saved
NAME= "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/meth.168.merged.20miss.txt"
N = 168
MIN.COV = 20
MAX.COV= 99.9


# Metadata for all samples
info <- read_xlsx("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/NCBI_info_updated_20220727.xlsx")

directory="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_output_files_test"

filenames <- list.files(path = directory, 
                        pattern = "5mC.cov\\.gz$", 
                        full.names = TRUE)
filenames=as.list(filenames)
names=list.files(path=directory,
                        pattern = "5mC.cov\\.gz$", 
                        full.names = TRUE)
names=gsub("__",".", names)
names=gsub("\\..*", "", names)
names=as.list(names)

my.methRaw=methRead(location = filenames,
                    sample.id = names,
                    assembly = "hirundo.rustica",
                    pipeline = 'bismarkCoverage',
                    context = "CpG",
                    treatment = c(rep(1,N)), # This is a fake treatment lable, no real meaning at all
                    mincov = 10)

filtered.my.methRaw <- filterByCoverage(my.methRaw, 
                                        lo.count=MIN.COV,
                                        lo.perc = NULL,
                                        hi.count = NULL,
                                        hi.perc = MAX.COV)

# normalize read coverage's between samples to avoid bias introduced by systematically more sequenced samples
normalized.myobj=normalizeCoverage(filtered.my.methRaw, method="median")

# merging samples for sites that are shared between them
meth=unite(normalized.myobj, min.per.group = 135L, destrand = F)

##number will change each time depending on the number of sites with 0 missing data
chr.names <- as.data.frame(unique(meth$chr))
meth_nosex <- meth[meth$chr != "chrW" & meth$chr != "chrW_random1" & meth$chr != "chrW_random2" & meth$chr != "chrZ" & meth$chr != "chrZ_random1", ] 

#get rid of 'random' chromosomes
meth_nosex_norand <- meth_nosex[meth_nosex$chr != "chr01A_random1" & 
                                  meth_nosex$chr != "chr01_random1" & 
                                  meth_nosex$chr != "chr02_random1" & 
                                  meth_nosex$chr != "chr02_random2" & 
                                  meth_nosex$chr != "chr02_random3" & 
                                  meth_nosex$chr != "chr02_random4" & 
                                  meth_nosex$chr != "chr03_random1" & 
                                  meth_nosex$chr != "chr03_random2" & 
                                  meth_nosex$chr != "chr03_random3" & 
                                  meth_nosex$chr != "chr03_random4" & 
                                  meth_nosex$chr != "chr03_random5" & 
                                  meth_nosex$chr != "chr04A_random1"& 
                                  meth_nosex$chr != "chr04_random1" & 
                                  meth_nosex$chr != "chr04_random2" & 
                                  meth_nosex$chr != "chr04_random3" & 
                                  meth_nosex$chr != "chr04_random4" & 
                                  meth_nosex$chr != "chr05_random1" & 
                                  meth_nosex$chr != "chr10_random1" & 
                                  meth_nosex$chr != "chr11_random1" & 
                                  meth_nosex$chr != "chr14_random1" & 
                                  meth_nosex$chr != "chr21_random1" & 
                                  meth_nosex$chr != "chr25_random1" & 
                                  meth_nosex$chr != "chr25_random2" & 
                                  meth_nosex$chr != "chr27_random1" & 
                                  meth_nosex$chr != "chr27_random2" & 
                                  meth_nosex$chr != "chr29_random1" & 
                                  meth_nosex$chr != "chr33_random1", ] 
meth.final <- meth_nosex_norand

write.table(meth.final, file = NAME, quote = FALSE, sep = "\t",
            row.names = FALSE)

