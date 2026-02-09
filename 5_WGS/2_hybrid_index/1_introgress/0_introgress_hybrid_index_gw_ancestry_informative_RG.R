#Script to run introgress on genome-wide data from barn swallow individuals 
# --- 1. Load necessary libraries ---
# Ensure these libraries are installed with install.packages() if you haven't already.
library(vcfR, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
library(introgress, lib = "/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2")
library(dplyr)
library(stringr)
library(adegenet)
library(ggplot2)
library(memuse)

# --- 2. Load VCF and Sample Metadata ---
# Load the main VCF file containing all samples.
vcf <- read.vcfR("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/2_hybrid_index/0_ancestry_informative/RGvcf.ancestry.informative.fst8.vcf.gz")

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

#Put this into a list
vcf_subsets <- list(all_vcf = all_vcf)

#-------5. Initial Data Preparation ---------
# Extract the genotype data of the ancestry informative SNPs
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
geno_subsets_numeric$all_vcf[1:5,1:5]

geno_subsets$all_vcf[1:5,1:5]

#----------6. calculate allele frequencies of parental populations------------
# Define parental sample index ranges (adjust if needed!)
G_idx <- c(1:5,7,9:12)     # columns 1–12 = G parental samples
R_idx <- c(60:63,65:67,69,70)    # columns 59–70 = R parental samples

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


# --- 7. Identify fixed SNPs and subset matrices per site.id ---
# ---- This part is redundant in this script because we have already chosen ancestry informative SNPs, so they already show differences in allele freq -----
find_fixed_snps <- function(geno_mat, af_df, cutoff = 0.4) {
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

# -----------11. Run introgress to determine HI, interspecific heterozygosity and parental allel frequencies ------

# Pull data for this site (for all individuals)
gen.mat <- fixed_subsets$all_vcf$gen.mat[, 1:70]
af_df   <- fixed_subsets$all_vcf$conv.mat
locus.info <- locus_data$all_vcf

# Skip if no SNPs or only 1 SNP
if (is.null(gen.mat) || nrow(gen.mat) <= 1) return(NULL)

# Get parental allele frequencies (here I use the raw matricies)
G.af <- fixed_subsets$all_vcf$gen.mat[,c(1:5,7,9:12)]
R.af <- fixed_subsets$all_vcf$gen.mat[,c(60:63,65:67,69,70)]

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

het.sim <- introgress::calc.intersp.het(
  introgress.data = count.matrix
)

hi.results <- introgress::est.h(
  introgress.data = count.matrix,
  loci.data       = locus.info
)

#double check heterozygosity calculation
het.manual <- colMeans(gen.mat == "0/1", na.rm = TRUE)

hi.het.hyb <- hi.results
hi.het.hyb$het <- het.sim
hi.het.hyb$het.manual <- het.manual
hi.het.hyb$sample <- all_samples

# ----------------Addition from Reviewers to look at triangle plots of HI and interspecific heterozygosity -----------------

#Run triangle plot from introgress
introgress::triangle.plot(hi.index=hi.results, int.het=het.sim, pdf=FALSE)

#Make nicer plotting in ggplot
plot_data <- hi.het.hyb %>%
  dplyr::rename(h_index = h, heterozygosity = het) %>% # Using dplyr::rename as requested
  mutate(color_group = case_when(
    grepl("_RG_", sample) ~ "RG",  # Check for RG first so it isn't caught by R or G alone
    grepl("_G_",  sample) ~ "G",
    grepl("_R_",  sample) ~ "R",
    TRUE                  ~ "Other"
  ))

# 2. Plot with manual color scales
ggplot(plot_data, aes(x = h_index, y = heterozygosity, color = color_group)) +
  geom_point(alpha = 0.7, size = 2) +
  # These segments represent the theoretical maximum heterozygosity
  geom_segment(x = 0, y = 0, xend = 0.5, yend = 1, linetype = "dashed", color = "black") +
  geom_segment(x = 1, y = 0, xend = 0.5, yend = 1, linetype = "dashed", color = "black") +
  # Apply your specific color panel
  scale_color_manual(values = c(
    "R"  = "#aa0a0f", 
    "G"  = "#6294b7", 
    "RG" = "#9716a8"
  )) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Hybrid index",
    y = "Interspecific heterozygosity",
    color = "Sample Type"
  )

# -----------12. Clean results and export for downstream use ------
write.table(hi.het.hyb, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowWGS_168/2_hybrid_index/1_introgress/1_RG_ancestry_informative_HI.txt", quote=F)



