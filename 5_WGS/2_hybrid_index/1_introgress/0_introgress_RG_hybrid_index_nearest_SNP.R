## ----- Use Introgress to get nearest SNP for Hybrid Models -------
## ---- Sarah Mueller - DNA methylation in Barn Swallows

# --- 1. Load necessary libraries ---
# Ensure these libraries are installed with install.packages() if you haven't already.
library(vcfR, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
library(introgress, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
library(dplyr)
library(stringr)
library(adegenet)
library(ggplot2)

# --- 2. Load VCF and Sample Metadata ---
# Load the main VCF file containing all samples.
vcf <- read.vcfR("/dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/HR.168.biallelic.maf05.dp_3_60.vcf.gz")

# Load the table containing sample names and population information.
sample_names_rg <- read.table("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.70.RGhyb.sub.pop.txt", header=T)

# Define the full list of samples to keep from the metadata file.
all_samples <- sample_names_rg$V1

# Define the hybrid sample names by filtering the metadata.
hybrid_samples <- sample_names_rg %>%
  filter(str_detect(V1, "_RG_")) %>%
  pull(V1) # Use pull() to extract the vector of names.

# --- 3. Subset the VCF files ---
# Get the sample names from the VCF's column headers.
# The first column, 'FORMAT', should be excluded.
vcf_sample_names <- colnames(vcf@gt)[-1]

# Create a VCF object containing all samples of interest.
all_vcf <- vcf[, c(TRUE, vcf_sample_names %in% all_samples)]

##FILTER FOR MAF 0.1 --to avoid monomorphic SNPS
maf <- maf(all_vcf, element=2)
maf.loci <- maf[,4]
keep_snps <- which(maf.loci >= 0.1)
all_vcf_filtered <- all_vcf[keep_snps, ]

#----------4. Load in candidate CpG regions/sites to get local ancestry from-----------
cpg_sites <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/4_RGT_Hybrids/3_methylation_models/RG_hyb/RG_diff/0_RG_parental_hybrid_tukey_results.txt", header=T)

cpg_sites <- cpg_sites %>%
  separate(site.id, into = c( "chr", "start", "end"), sep = "_", remove = FALSE)

cpg_sites <- cpg_sites %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  )

#-------5. Loop through unique site.ids and make a list with individual vcfs for the nearest SNP-------------
# Create an empty list to store VCF subsets
vcf_subsets <- list()

# SNP positions and chromosomes from VCF
snp_pos <- as.numeric(all_vcf_filtered@fix[,2])
snp_chr <- all_vcf_filtered@fix[,1]

# Loop through each CpG site
for (i in 1:nrow(cpg_sites)) {
  
  site_id <- cpg_sites$site.id[i]
  chr     <- cpg_sites$chr[i]
  start   <- cpg_sites$start[i]
  end     <- cpg_sites$end[i]
  
  # If end is NA, use start as CpG position
  # Otherwise, use midpoint of region
  if (is.na(end)) {
    cpg_pos <- start
  } else {
    cpg_pos <- round((start + end) / 2)
  }
  
  # Subset SNPs on the same chromosome
  snp_chr_idx <- which(snp_chr == chr)
  
  if (length(snp_chr_idx) > 0) {
    snp_chr_pos <- snp_pos[snp_chr_idx]
    
    # Find nearest SNP by absolute distance
    nearest_idx <- snp_chr_idx[which.min(abs(snp_chr_pos - cpg_pos))]
    
    # Store the nearest SNP VCF row
    vcf_subsets[[site_id]] <- all_vcf_filtered[nearest_idx, ]
  } else {
    message("No SNPs found on chromosome ", chr, " for site ", site_id)
  }
}



#---------6.Extract genotxype information from nearest SNP loci------------
# Convert the vcf_subsets list into a tidy dataframe
geno_table <- lapply(names(vcf_subsets), function(site_id) {
  vcf <- vcf_subsets[[site_id]]
  
  # Extract chromosome and position
  chr <- vcf@fix[1, "CHROM"]
  pos <- as.numeric(vcf@fix[1, "POS"])
  
  # Extract genotype calls (matrix: sites x individuals)
  gt <- extract.gt(vcf)  # gives matrix, rows=sites (1 here), cols=individuals
  gt_vec <- as.vector(gt[1, ])  # get first row (only one SNP)
  names(gt_vec) <- colnames(gt) # individual IDs as names
  
  # Build a single row as dataframe
  data.frame(
    site_id = site_id,
    chr = chr,
    pos = pos,
    t(gt_vec),        # transpose vector → individuals become columns
    check.names = FALSE
  )
})

# Bind all rows together
geno_table <- bind_rows(geno_table)



# Function to convert GT strings into numeric dosage
gt_to_dosage <- function(gt) {
  gt <- gsub("\\|", "/", gt)  # make sure phased (|) → unphased (/)
  sapply(gt, function(x) {
    if (x %in% c("0/0")) return(0)
    else if (x %in% c("0/1", "1/0")) return(1)
    else if (x %in% c("1/1")) return(2)
    else return(NA)  # for "./." or missing
  })
}

geno_table <- lapply(names(vcf_subsets), function(site_id) {
  vcf <- vcf_subsets[[site_id]]
  
  # SNP info
  chr <- vcf@fix[1, "CHROM"]
  pos <- as.numeric(vcf@fix[1, "POS"])
  
  # Extract genotypes
  gt <- extract.gt(vcf)
  gt_vec <- gt_to_dosage(gt[1, ])  # convert to dosage
  
  # Build row
  data.frame(
    site_id = site_id,
    chr = chr,
    pos = pos,
    t(gt_vec),
    check.names = FALSE
  )
})

# Final combined table
geno_table <- bind_rows(geno_table)


write.table(geno_table, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/1_introgress/1_RG_local_nearest_SNP_maf1.txt", row.names=F, quote=F)
