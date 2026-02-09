## ----- Use Introgress to get local (cis-) ancestry blocks for Hybrid Models -------
## ---- Sarah Mueller - DNA methylation in Barn Swallows

# --- 1. Load necessary libraries ---
# Ensure these libraries are installed with install.packages() if you haven't already.
library(vcfR)
library(introgress)
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

#----------4. Load in candidate CpG regions/sites to get local ancestry from-----------
cpg_sites <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/4_RGT_Hybrids/3_methylation_models/RG_hyb/RG_diff/0_RG_parental_hybrid_tukey_results.txt", header=T)

cpg_sites <- cpg_sites %>%
  separate(site.id, into = c( "chr", "start", "end"), sep = "_", remove = FALSE)

cpg_sites <- cpg_sites %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  )

#-------5. Loop through unique site.ids and make a list with individual vcfs for each local region-------------
# Create an empty list to store VCF subsets
vcf_subsets <- list()

# Loop through each CpG site
for (i in 1:nrow(cpg_sites)) {
  
  site_id <- cpg_sites$site.id[i]
  chr <- cpg_sites$chr[i]
  start <- cpg_sites$start[i]
  end <- cpg_sites$end[i]
  
  # If end is NA, define region as start ± 40kb
  if (is.na(end)) {
    region_start <- max(1, start - 5000)   # avoid negative coordinates
    region_end   <- start + 5000
  } else {
    region_start <- max(1, start - 5000)
    region_end   <- end + 5000
  }
  
  # Extract SNPs from VCF in this region
  # Note: vcfR::extract.gt uses chrom/pos stored in @fix slot
  snp_pos <- as.numeric(vcf@fix[,2])   # positions
  snp_chr <- all_vcf@fix[,1]               # chromosomes
  
  idx <- which(snp_chr == chr & snp_pos >= region_start & snp_pos <= region_end)
  
  if (length(idx) > 0) {
    vcf_region <- all_vcf[idx, ]
    vcf_subsets[[site_id]] <- vcf_region
  } else {
    message("No SNPs found for ", site_id, " (", chr, ":", region_start, "-", region_end, ")")
  }
}

# --------- 6. Initial Data Preparation ---------
# Extract the genotype data for each local ancestry dataframe
geno_subsets <- lapply(vcf_subsets, function(vcf_region) {
  extract.gt(vcf_region, element = "GT", as.numeric = FALSE)
})

#Convert genotype strings to numeric codes
geno_subsets_numeric <- lapply(geno_subsets, function(conv.mat) {
  # Replace genotype strings with numeric codes
  conv.mat[conv.mat == "0/0"] <- 0
  conv.mat[conv.mat == "0/1"] <- 1
  conv.mat[conv.mat == "1/1"] <- 2
  
  # Convert to numeric while keeping matrix structure
  conv.mat <- matrix(as.numeric(conv.mat), 
                     nrow = nrow(conv.mat), 
                     ncol = ncol(conv.mat),
                     dimnames = dimnames(conv.mat))
  
  return(conv.mat)
})

##Test a few outputs
geno_subsets_numeric$chr01_105650270_105652270[1:5,1:5]

#----------7. calculate allele frequencies of parental populations------------
# Define parental sample index ranges (adjust if needed!)
G_idx <- 1:12     # columns 1–12 = G parental samples
R_idx <- 59:70    # columns 59–70 = R parental samples

# Function to calculate AF given a genotype matrix
calc_AF <- function(geno_mat) {
  # G parental AF
  G.af <- (rowSums(geno_mat[, G_idx, drop=FALSE], na.rm = TRUE) /
             rowSums(!is.na(geno_mat[, G_idx, drop=FALSE]))) / 2
  
  # R parental AF
  R.af <- (rowSums(geno_mat[, R_idx, drop=FALSE], na.rm = TRUE) /
             rowSums(!is.na(geno_mat[, R_idx, drop=FALSE]))) / 2
  
  # Return a dataframe with both
  data.frame(G.af = G.af, R.af = R.af)
}

# Apply to each entry in geno_subsets
af_subsets <- lapply(geno_subsets_numeric, calc_AF)


# --- 8. Identify fixed SNPs and subset matrices per site.id ---
find_fixed_snps <- function(geno_mat, af_df, cutoff = 0.1) {
  # Calculate diff between parental AF
  diff <- abs(af_df$G.af - af_df$R.af)
  
  # Logical index of fixed SNPs
  fixed_idx <- !is.na(diff) & diff > cutoff
  
  # Count fixed vs non-fixed
  counts <- table(fixed_idx)
  
  # Subset genotype matrix (gen.mat)
  gen.mat <- geno_mat[fixed_idx, , drop = FALSE]
  
  # Subset AF dataframe (conv.mat equivalent)
  conv.mat <- af_df[fixed_idx, , drop = FALSE]
  
  return(list(
    diff = diff,
    counts = counts,
    gen.mat = gen.mat,
    conv.mat = conv.mat
  ))
}

# Apply to all site.id entries
fixed_subsets <- mapply(
  FUN = find_fixed_snps,
  geno_mat = geno_subsets,
  af_df = af_subsets,
  SIMPLIFY = FALSE
)

# --- 9. Find site.ids with at least one fixed SNP ---
sites_with_fixed <- names(fixed_subsets)[
  sapply(fixed_subsets, function(x) {
    "TRUE" %in% names(x$counts) && x$counts[["TRUE"]] > 0
  })
]

# Print them
sites_with_fixed

# --- 10. Build locus_data list from fixed_subsets ---
locus_data <- lapply(fixed_subsets, function(x) {
  gen.mat <- x$gen.mat
  
  if (nrow(gen.mat) > 0) {
    locus.info <- data.frame(
      locus = rownames(gen.mat),
      type  = rep("C", times = nrow(gen.mat))
    )
  } else {
    locus.info <- data.frame(
      locus = character(0),
      type  = character(0)
    )
  }
  
  locus.info
})

# Keep the site.id names
names(locus_data) <- names(fixed_subsets)

# --- 11. Run introgress hybrid index estimation per site.id ---
run_introgress_fixed <- function(site_id, fixed_subsets, locus_data) {
  # Pull data for this site
  gen.mat <- fixed_subsets[[site_id]]$gen.mat
  af_df   <- fixed_subsets[[site_id]]$conv.mat
  locus.info <- locus_data[[site_id]]
  
  # Skip if no SNPs or only 1 SNP
  if (is.null(gen.mat) || nrow(gen.mat) <= 1) return(NULL)
  
  # Get parental allele frequencies
  G.af <- af_df$G.af
  R.af <- af_df$R.af
  
  # Prepare introgress data
  count.matrix <- introgress::prepare.data(
    admix.gen  = gen.mat,
    loci.data  = locus.info,
    parental1  = G.af,
    parental2  = R.af,
    pop.id     = FALSE,
    ind.id     = FALSE,
    fixed      = FALSE
  )
  
  # Estimate hybrid index
  hi.index.sim <- introgress::est.h(
    introgress.data = count.matrix,
    loci.data       = locus.info
  )
  
  return(hi.index.sim)
}

# Apply to all site_ids
hi_results_fixed <- lapply(names(fixed_subsets), run_introgress_fixed,
                           fixed_subsets = fixed_subsets,
                           locus_data    = locus_data)

# Keep names consistent
names(hi_results_fixed) <- names(fixed_subsets)


# -----------12. Clean results and put together to run models ------
# Extract h columns and rename to site IDs
h_list <- lapply(names(hi_results_fixed), function(site_id) {
  hi_df <- hi_results_fixed[[site_id]]
  
  if(!"h" %in% colnames(hi_df)) return(NULL)  # skip if no h
  
  h_col <- hi_df[,"h", drop = FALSE]          # keep as data.frame
  colnames(h_col) <- site_id                  # rename column to site ID
  h_col
})

# Remove any NULLs
h_list <- h_list[!sapply(h_list, is.null)]

# Combine all columns into a single data.frame
h_combined <- do.call(cbind, h_list)

h_combined <- as.data.frame(t(h_combined))
colnames(h_combined) <- colnames(geno_subsets$chr03_10498_11502)

write.table(h_combined, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/1_introgress/1_RG_local_ancestry_HI.txt", quote=F)
