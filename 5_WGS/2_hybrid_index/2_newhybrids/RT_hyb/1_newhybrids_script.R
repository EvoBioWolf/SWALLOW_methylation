###Newhybrid runs

# --- 1. Load necessary libraries ---
# Ensure these libraries are installed with install.packages() if you haven't already.
.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))
library(vcfR, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
library(introgress, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
library(dplyr)
library(stringr)
library(adegenet)
library(ggplot2)
library(hybriddetective)
library(genepopedit)
library(parallelnewhybrid)

# --- 2. Load VCF and Sample Metadata ---
# Load the main VCF file containing all samples.
vcf <- read.vcfR("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/0_ancestry_informative/RT.ancestry.informative.fst8.vcf.gz")

# Load the table containing sample names and population information.
sample_names_rg <- read.table("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.60.RThyb.sub.pop.txt", header=T)

# Define the full list of samples to keep from the metadata file.
all_samples <- sample_names_rg$V1

# Define the hybrid sample names by filtering the metadata.
hybrid_samples <- sample_names_rg %>%
  filter(str_detect(V1, "_RG_")) %>%
  pull(V1) # Use pull() to extract the vector of names.

# The first column, 'FORMAT', should be excluded.
vcf_sample_names <- colnames(vcf@gt)[-1]

# Create a VCF object containing all samples of interest.
all_vcf <- vcf[, c(TRUE, vcf_sample_names %in% all_samples)]

# -----3. Extract information from the VCF file ----- 
gt <- extract.gt(all_vcf, element = "GT", as.numeric = FALSE)

# -----4. Change the allele formatting to fit the new hybrids software --- 
allele_map <- function(gt) {
  if (is.na(gt) || gt %in% c("./.", ".|.")) return("000000")
  if (gt %in% c("0/0","0|0")) return("100100")
  if (gt %in% c("0/1","1/0","0|1","1|0")) return("100110")
  if (gt %in% c("1/1","1|1")) return("110110")
  return("000000")
}

nh_mat <- apply(gt, c(1,2), allele_map)
nh_mat <- t(nh_mat)   # change to ind = rows, loci = columns

# ----- 5. Add correct heading ----
num_ind  <- nrow(nh_mat)
num_loci <- ncol(nh_mat)

header <- c(
  paste("NumIndivs", num_ind),
  paste("NumLoci", num_loci),
  "Digits 3",
  "Format Lumped",
  paste("LocusNames", paste(rownames(gt), collapse = " "))
)


# ---- 6. Write input file ------------
out <- file("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RT_hyb/2_empirical/RT_60_NH.txt", "w")

writeLines(header, out)

for (i in seq_len(num_ind)) {
  line <- paste(i, paste(nh_mat[i, ], collapse = " "))
  writeLines(line, out)
}

close(out)

# ----- 6.5 Write sample names input file
df_inds <- data.frame(ID = colnames(gt))

write.table(
  df_inds,
  file = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RT_hyb/2_empirical/RT_60_NH_individuals.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# ----- 7. Subset columns to create ref_set -----
# Only include individuals who were identified are pure parentals from introgress
# For this check the h values in ./5_RGT_Hybrids/3_methylation_boxplots/1_RG_ancestry_informative_HI.txt
# This forms the basis of simulations so we can tell how good the loci are at defining hybrid classes

#ref_set <- nh_mat[c(1:4,7:12,50:54,56,59,60),]

#write_genepop <- function(mat, file) {
#  
#  loci <- colnames(mat)
#  inds <- rownames(mat)
  
#  popG <- inds[1:10]
#  popR <- inds[11:length(inds)]
  
#  con <- file(file, "w")
  
#  writeLines("Reference parentals", con)
#  writeLines(loci, con)
  
#  writeLines("POP", con)
#  for (i in popG) {
#    writeLines(paste(i, " ,  ", paste(mat[i, ], collapse = " "), sep = ""), con)
#  }
  
#  writeLines("POP", con)
#  for (i in popR) {
#    writeLines(paste(i, " ,  ", paste(mat[i, ], collapse = " "), sep = ""), con)
#  }
  
#  close(con)
#}

#write_genepop(
#  mat  = ref_set,
#  file = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RT_hyb/1_simulations/RT_par_sim.txt"
#)


# ---- 7. From the parental dataset, generate simulated data of hybrids to see how good loci are at detecting hybrids
## This will create three independent simulations, each replicated three times.
#hybriddetective::freqbasedsim_GTFreq(GenePopData = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RT_hyb/1_simulations/RT_par_sim.txt", pop.groups = c("PopR", "PopT"), sample.sizePure = 11, sample.sizeF1 = 10, sample.sizeF2 = 10, sample.sizeBC = 10, NumSims = 3, NumReps = 3)


# ----- 8. Run NewHybrids on empirical and simulated data -------
# file path to the folder in which NewHybrids is installed. Note: this folder must be named "newhybrids"
your.NH <- "/dss/dsslegfs01/pr53da/pr53da-dss-0034/newhybrids"

# Execute parallelnh. NOTE: "xx" must be replaced by the correct designation for your operating system.
# Here I run on the simulated data
parallelnewhybrid::parallelnh_LINUX(folder.data = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RT_hyb/1_simulations/", where.NH = your.NH, burnin = 100000, sweeps = 200000)

#Run again on empirical data with same parameters
#parallelnewhybrid::parallelnh_LINUX(folder.data = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RT_hyb/2_empirical/", where.NH = your.NH, burnin = 100000, sweeps = 200000)

