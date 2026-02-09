## LD decay plotting script
library(tidyverse)
library(fuzzyjoin)

#Examine if theres any overlap in the SNP and methylation files
snp.meth <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_Rall/HR_Rall_snp.meth.maf2.ldwin100.gz", header= T)

snp.meth.ld <- snp.meth %>%
  mutate(distance=abs(BP_A - BP_B))

snp.meth.ld <- snp.meth.ld %>%
  filter(!distance == 0) %>%
  filter(!(distance == 1 & startsWith(SNP_A, "M_")))

snp.meth.ld.70kb <- snp.meth.ld %>%
  mutate(distance_bin = cut(distance, 
                            breaks = seq(0, 70000, by = 100), 
                            labels = seq(100, 70000, by = 100), 
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%  # Convert factor to numeric
  group_by(distance_bin) %>%
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

snp.meth.ld.5kb <- snp.meth.ld %>%
  mutate(distance_bin = cut(distance, 
                            breaks = seq(0, 5000, by = 20), 
                            labels = seq(20, 5000, by = 20), 
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%  # Convert factor to numeric
  group_by(distance_bin) %>%
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

ggplot(snp.meth.ld.70kb, aes(x = distance_bin, y = avg_R2)) +
  geom_point(alpha = 0.5, color = "orange") +  # Scatterplot of points
  geom_smooth(method = "loess", se = TRUE, color = "darkorange", span = 0.2) +  # Smoothed trend line
  #facet_wrap(~chr)+ 
  labs(title = "SNP-METH LD Decay Plot",
       x = "Distance",
       y = "Average R²") +
  theme_minimal()

ggplot(snp.meth.ld.5kb, aes(x = distance_bin, y = avg_R2), na.rm=F) +
  geom_point(alpha = 0.5, color = "orange") +  # Scatterplot of points
  #geom_smooth(method = "loess", se = TRUE, color = "darkorange", span = 0.2) +  # Smoothed trend line
  #facet_wrap(~chr)+ 
  labs(title = "SNP-CpG LD Decay Plot",
       x = "Distance",
       y = "Average R²") +
  theme_minimal()


##Examine this first bin of SNP-Meth LD
snp.meth.dist20 <- snp.meth.ld %>%
  filter(distance < 20)

##Write these sites out to a table, pull them from the VCF and see the end to look in-depth at these SNP-Meth comparisons
#write.table(dist1, "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_Rall/HR.85.dist1.txt", quote=F, row.names=F)

#ONLY SNP LD
snp.ld <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_Rall/HR_Rall_snp.maf2.ldwin100.gz", header= F)
colnames(snp.ld) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B",  "BP_B", "SNP_B", "R2") #colnames(snp.meth.ld[1:7])

snp.ld <- snp.ld %>%
  mutate(distance=abs(BP_A - BP_B))

snp.ld.70kb <- snp.ld %>%
  mutate(distance_bin = cut(distance, 
                            breaks = seq(0, 70000, by = 100), 
                            labels = seq(100, 70000, by = 100), 
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%  # Convert factor to numeric
  group_by(distance_bin) %>%
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

snp.ld.5kb <- snp.ld %>%
  mutate(distance_bin = cut(distance, 
                            breaks = seq(0, 5000, by = 20), 
                            labels = seq(20, 5000, by = 20), 
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%  # Convert factor to numeric
  group_by(distance_bin) %>%
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

ggplot(snp.ld.70kb, aes(x = distance_bin, y = avg_R2)) +
  geom_point(alpha = 0.5, color = "blue") +  # Scatterplot of points
  #geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.2) +  # Smoothed trend line
  #facet_wrap(~chr)+ 
  labs(title = "SNP LD Decay Plot",
       x = "Distance",
       y = "Average R²") +
  theme_minimal()

ggplot(snp.ld.5kb, aes(x = distance_bin, y = avg_R2), na.rm=F) +
  geom_point(alpha = 0.5, color = "grey22") +  # Scatterplot of points
  #geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.2) +  # Smoothed trend line
  #facet_wrap(~chr)+ 
  labs(title = "SNP LD Decay Plot",
       x = "Distance",
       y = "Average R²") +
  theme_minimal()

###ONLY METHYLATION LD
meth.ld <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_Rall/HR_Rall_meth.maf2.ldwin100.gz", header= F)
colnames(meth.ld) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B",  "BP_B", "SNP_B", "R2")

meth.ld <- meth.ld %>%
  mutate(distance=abs(BP_A - BP_B))

meth.ld.70kb <- meth.ld %>%
  mutate(distance_bin = cut(distance, 
                            breaks = seq(0, 70000, by = 100), 
                            labels = seq(100, 70000, by = 100), 
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%  # Convert factor to numeric
  group_by(distance_bin) %>%
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

meth.ld.5kb <- meth.ld %>%
  mutate(distance_bin = cut(distance, 
                            breaks = seq(0, 5000, by = 20), 
                            labels = seq(20, 5000, by = 20), 
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%  # Convert factor to numeric
  group_by(distance_bin) %>%
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

ggplot(meth.ld.70kb, aes(x = distance_bin, y = avg_R2)) +
  geom_point(alpha = 0.5, color = "blue") +  # Scatterplot of points
  #geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.2) +  # Smoothed trend line
  #facet_wrap(~chr)+ 
  labs(title = "METH LD Decay Plot",
       x = "Distance",
       y = "Average R²") +
  theme_minimal()

ggplot(meth.ld.5kb, aes(x = distance_bin, y = avg_R2), na.rm=F) +
  geom_point(alpha = 0.5, color = "blue") +  # Scatterplot of points
  #geom_smooth(method = "loess", se = FALSE, color = "red", span = 0.2) +  # Smoothed trend line
  #facet_wrap(~chr)+ 
  labs(title = "CpG LD Decay Plot",
       x = "Distance",
       y = "Average R²") +
  theme_minimal()

####Combine into one graph

df_list <- list(snp.meth = snp.meth.ld.5kb, snp = snp.ld.5kb, meth = meth.ld.5kb)

snp.meth.combo.5kb <- bind_rows(lapply(names(df_list), function(name) {
  df_list[[name]] %>% mutate(Type = name)
}))

print(snp.meth.combo.5kb)

ggplot(snp.meth.combo.5kb, aes(x = distance_bin, y = avg_R2, color= Type), na.rm=F) +
  geom_point(alpha = 0.5, size=2) +  # Scatterplot of points
  #geom_smooth(method = "loess", se = FALSE, span = 0.2) +  # Smoothed trend line
  scale_color_manual(values= c("#e53207", "grey33", "#48497e"))+
  labs(title = "Linkage (5kb)",
       x = "Distance",
       y = "Average R²") +
  theme_minimal() +
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_text(face = "bold", size=16),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_blank())


##Methylation LD plots
meth.info <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/1_CpG_islands/168.20miss/homer.loci.prop.168.merged.20miss.rmlowv.rmout.txt", header=T)
meth.df <- meth.info %>%
        separate(site.id, into = c("chr", "pos"), sep= "_")
meth.df$pos <- as.numeric(meth.df$pos)

##add region information and PST
pst_homer <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/1_PST/0_site_pst/0_pst_variant_c05.txt", header=T)
pst_homer_df <- pst_homer %>%
  separate(site.id, into = c("chr", "pos"), sep= "_")
pst_homer_df$pos <- as.numeric(pst_homer_df$pos)

meth.ld <- meth.ld %>%
  mutate(
    CHR_A = case_when(
      CHR_A %in% c('chr01A', 'chr04A') ~ CHR_A,
      CHR_A %in% as.character(1:9) ~ paste0('chr0', CHR_A),
      TRUE ~ paste0('chr', CHR_A)
    ),
    CHR_B = case_when(
      CHR_B %in% c('chr01A', 'chr04A') ~ CHR_B,
      CHR_B %in% as.character(1:9) ~ paste0('chr0', CHR_B),
      TRUE ~ paste0('chr', CHR_B)
    ),
    distance = (BP_B - BP_A)
  )

meth.merged.ld.a <- meth.ld %>%
      left_join(pst_homer_df, by = c("CHR_A" = "chr", "BP_A" = "pos"))
meth.merged.ld.a.b <- meth.merged.ld.a %>%
  left_join(pst_homer_df, by = c("CHR_B" = "chr", "BP_B" = "pos"))

#get just high ld sites
meth.merged.ld.a.b.highld <- meth.merged.ld.a.b %>%
      filter(R2 > 0.5)

# Count occurrences of total data
meth_ld_counts <- meth.merged.ld.a.b %>%
  dplyr::count(category.y)

#Count occurances in high LD data
meth_highld_counts <- meth.merged.ld.a.b.highld %>%
  dplyr::count(category.y)

# Pie Chart for `sea`
p2 <- ggplot(meth_ld_counts, aes(x = "", y = n, fill = category.y)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()  +
  theme(legend.position = "none")

p3 <- ggplot(meth_highld_counts, aes(x = "", y = n, fill = category.y)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()


p2 +p3

###JUST region and annotation
meth.merged.ld.a.b <- meth.merged.ld.a.b %>%
  separate(Annotation.x, into = c("sea", "region"), sep = "_", remove = FALSE)
meth.merged.ld.a.b.highld <- meth.merged.ld.a.b.highld %>%
  separate(Annotation.x, into = c("sea", "region"), sep = "_", remove = FALSE)

# Count occurrences of total data
meth_ld_sea <- meth.merged.ld.a.b %>%
  dplyr::count(sea)

#Count occurances in high LD data
meth_highld_sea <- meth.merged.ld.a.b.highld %>%
  dplyr::count(sea)

# Pie Chart for `sea`
p4 <- ggplot(meth_ld_sea, aes(x = "", y = n, fill = sea)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()  +
  theme(legend.position = "none")

p5 <- ggplot(meth_highld_sea, aes(x = "", y = n, fill = sea)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()


p4 +p5

# Count occurrences of total data
meth_ld_region <- meth.merged.ld.a.b %>%
  dplyr::count(region)

#Count occurances in high LD data
meth_highld_region <- meth.merged.ld.a.b.highld %>%
  dplyr::count(region)

# Pie Chart for `sea`
p6 <- ggplot(meth_ld_region, aes(x = "", y = n, fill = region)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()  +
  theme(legend.position = "none")

p7 <- ggplot(meth_highld_region, aes(x = "", y = n, fill = region)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()


p6 +p7


###Calculate bins based on region
meth.ld.5kb.category <- meth.merged.ld.a.b %>%
  mutate(distance_bin = cut(distance,
                            breaks = seq(0, 5000, by = 20),
                            labels = seq(20, 5000, by = 20),
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%
  group_by(distance_bin, category.x) %>% # Group by both distance_bin AND category
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")


meth.ld.5kb.category <- meth.ld.5kb.category %>%
          filter(!category.x %in% c("island_TTS", "island_exon", "island_intron", "opensea_TTS", "opensea_exon", "opensea_intron", "shelf_TTS", "shelf_exon", "shelf_intron", "shore_TTS", "shore_exon", "shore_intron", "NA"))

ggplot(meth.ld.5kb.category, aes(x = distance_bin, y = avg_R2, color= category.x), na.rm=F) +
  #geom_point(alpha = 0.5, size=2) +  # Scatterplot of points
  geom_smooth(method = "loess", se = FALSE, span = 0.2) +  # Smoothed trend line
  #scale_color_manual(values= c("#e53207", "grey33", "#48497e"))+
  labs(title = "Linkage (5kb)",
       x = "Distance",
       y = "Average R²") +
  theme_minimal() +
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_blank())
#############################################################################
#####SAME FOR SNP-METH
#get it into a workable df....without so many sites
snp.meth.ld <- snp.meth.ld %>%
  filter(distance < 20000)

snp.meth.ld <- snp.meth.ld %>%
  mutate(
    CHR_A = case_when(
      CHR_A %in% c('chr01A', 'chr04A') ~ CHR_A,
      CHR_A %in% as.character(1:9) ~ paste0('chr0', CHR_A),
      TRUE ~ paste0('chr', CHR_A)
    ),
    CHR_B = case_when(
      CHR_B %in% c('chr01A', 'chr04A') ~ CHR_B,
      CHR_B %in% as.character(1:9) ~ paste0('chr0', CHR_B),
      TRUE ~ paste0('chr', CHR_B)
    ),
    distance = (BP_B - BP_A)
  )

snp.meth.merged.ld.a <- snp.meth.ld %>%
  inner_join(pst_homer_df, by = c("CHR_A" = "chr", "BP_A" = "pos"))
snp.meth.merged.ld.b <- snp.meth.ld %>%
  inner_join(pst_homer_df, by = c("CHR_B" = "chr", "BP_B" = "pos"))

snp.meth.merged.ld.a.b <- rbind(snp.meth.merged.ld.a, snp.meth.merged.ld.b)

#get just high ld sites
snp.meth.merged.ld.a.b.highld <- snp.meth.merged.ld.a.b %>%
  filter(R2 > 0.5)

# Count occurrences of total data
snp_meth_ld_counts <- snp.meth.merged.ld.a.b %>%
  dplyr::count(category)

#Count occurances in high LD data
snp_meth_highld_counts <- snp.meth.merged.ld.a.b.highld %>%
  dplyr::count(category)

# Pie Chart for `sea`
p2 <- ggplot(snp_meth_ld_counts, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Total CpG Regions") +
  theme_void()  +
  theme(legend.position = "none")

p3 <- ggplot(snp_meth_highld_counts, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()


p2 +p3

###JUST region and annotation
snp.meth.merged.ld.a.b <- snp.meth.merged.ld.a.b %>%
  separate(Annotation, into = c("sea", "region"), sep = "_", remove = FALSE)
snp.meth.merged.ld.a.b.highld <- snp.meth.merged.ld.a.b.highld %>%
  separate(Annotation, into = c("sea", "region"), sep = "_", remove = FALSE)

# Count occurrences of total data
snp_meth_ld_sea <- snp.meth.merged.ld.a.b %>%
  dplyr::count(sea)

#Count occurances in high LD data
snp_meth_highld_sea <- snp.meth.merged.ld.a.b.highld %>%
  dplyr::count(sea)

# Pie Chart for `sea`
p4 <- ggplot(snp_meth_ld_sea, aes(x = "", y = n, fill = sea)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Total CpG Regions") +
  theme_void()  +
  theme(legend.position = "none")

p5 <- ggplot(snp_meth_highld_sea, aes(x = "", y = n, fill = sea)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()


p4 +p5

# Count occurrences of total data
snp_meth_ld_region <- snp.meth.merged.ld.a.b %>%
  dplyr::count(region)

#Count occurances in high LD data
snp_meth_highld_region <- snp.meth.merged.ld.a.b.highld %>%
  dplyr::count(region)

# Pie Chart for `sea`
p6 <- ggplot(snp_meth_ld_region, aes(x = "", y = n, fill = region)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Total CpG Regions") +
  theme_void()  +
  theme(legend.position = "none")

p7 <- ggplot(snp_meth_highld_region, aes(x = "", y = n, fill = region)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD CpG Regions") +
  theme_void()


p6 +p7


###Calculate bins based on region
snp.meth.ld.5kb.category <- snp.meth.merged.ld.a.b %>%
  mutate(distance_bin = cut(distance,
                            breaks = seq(0, 5000, by = 20),
                            labels = seq(20, 5000, by = 20),
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%
  group_by(distance_bin, category) %>% # Group by both distance_bin AND category
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

snp.meth.ld.5kb.category <- snp.meth.ld.5kb.category %>%
  filter(!category %in% c("island_TTS", "island_exon", "island_intron", "opensea_TTS", "opensea_exon", "opensea_intron", "shelf_TTS", "shelf_exon", "shelf_intron", "shore_TTS", "shore_exon", "shore_intron", "NA"))

ggplot(snp.meth.ld.5kb.category, aes(x = distance_bin, y = avg_R2, color= category), na.rm=F) +
  #geom_point(alpha = 0.5, size=2) +  # Scatterplot of points
  geom_smooth(method = "loess", se = FALSE, span = 0.2) +  # Smoothed trend line
  #scale_color_manual(values= c("#e53207", "grey33", "#48497e"))+
  labs(title = "Linkage (5kb)",
       x = "Distance",
       y = "Average R²") +
  theme_minimal() +
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_blank())

###Calculate bins based on region
snp.meth.ld.5kb.sea <- snp.meth.merged.ld.a.b %>%
  mutate(distance_bin = cut(distance,
                            breaks = seq(0, 5000, by = 20),
                            labels = seq(20, 5000, by = 20),
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%
  group_by(distance_bin, sea) %>% # Group by both distance_bin AND category
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

ggplot(snp.meth.ld.5kb.sea, aes(x = distance_bin, y = avg_R2, color= sea), na.rm=F) +
  #geom_point(alpha = 0.5, size=2) +  # Scatterplot of points
  geom_smooth(method = "loess", se = FALSE, span = 0.2) +  # Smoothed trend line
  #scale_color_manual(values= c("#e53207", "grey33", "#48497e"))+
  labs(title = "Linkage (5kb)",
       x = "Distance",
       y = "Average R²") +
  theme_minimal() +
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_blank())

###Calculate bins based on region
snp.meth.ld.5kb.region <- snp.meth.merged.ld.a.b %>%
  mutate(distance_bin = cut(distance,
                            breaks = seq(0, 5000, by = 20),
                            labels = seq(20, 5000, by = 20),
                            right = TRUE)) %>%
  mutate(distance_bin = as.numeric(as.character(distance_bin))) %>%
  group_by(distance_bin, region) %>% # Group by both distance_bin AND category
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

ggplot(snp.meth.ld.5kb.region, aes(x = distance_bin, y = avg_R2, color= region), na.rm=F) +
  #geom_point(alpha = 0.5, size=2) +  # Scatterplot of points
  geom_smooth(method = "loess", se = FALSE, span = 0.2) +  # Smoothed trend line
  #scale_color_manual(values= c("#e53207", "grey33", "#48497e"))+
  labs(title = "Linkage (5kb)",
       x = "Distance",
       y = "Average R²") +
  theme_minimal() +
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_blank())

snp.meth.merged.ld.a.b.ld <- snp.meth.merged.ld.a.b %>%
  filter(R2 > 0.2)

ggplot(snp.meth.merged.ld.a.b.ld, aes(x = CHR_A, y = R2), na.rm=T) +
  geom_boxplot(coef = 1, outlier.shape = NA, fill= "skyblue")+
  labs(title = "Linkage by chromsome",
       y = "Average R²") +
  theme_minimal() +
  coord_flip()+
  ylim(0.2,0.7)+
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_blank(),
        plot.title = element_blank())


meth.merged.ld.a.b.ld <- meth.merged.ld.a.b %>%
  filter(R2 > 0.1)

ggplot(meth.merged.ld.a.b.ld, aes(x = CHR_A, y = R2), na.rm=T) +
  geom_boxplot(coef = 1, outlier.shape = NA, fill= "skyblue")+
  labs(title = "Linkage by chromsome",
       y = "Average R²") +
  theme_minimal() +
  coord_flip()+
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_blank(),
        plot.title = element_blank())

ggplot(chr11_data, aes(x = as.factor(BP_A), y = R2), na.rm=T) +
  geom_boxplot()+
  labs(title = "Linkage by chromsome",
       x = "Chromsome",
       y = "Average R²") +
  theme_minimal() +
  theme(legend.key.size = unit(0.75, "cm"), # Adjust size as needed
        legend.text = element_text(size = 14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(fill = "white", color = "black"), 
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_blank())


# Filter for chromosome 11
chr11_data <- snp.meth.merged.ld.a.b.ld %>%
  filter(CHR_A == "chr11" & CHR_B == "chr11") # Ensure both SNPs are on chr11

# Create the heatmap
ggplot(chr11_data, aes(x = as.factor(BP_A), y = as.factor(BP_B), fill = R2)) +
  geom_tile() +  # Use geom_tile for the heatmap cells
  scale_fill_gradient(low = "red", high = "blue") +  # Customize color scale
  labs(title = "Linkage Disequilibrium (R²) on Chromosome 11",
       x = "Base Position A",
       y = "Base Position B",
       fill = "R²") +
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        panel.grid = element_blank()) +  # Remove grid lines
  coord_equal() #Ensure that the x and y axis are equally scaled. this is important for a heatmap.

# 1. Create 50kb bins for BP_A and BP_B
chr11_data_binned <- snp.meth.merged.ld.a.b.ld %>%
  filter(CHR_A == "chr11" & CHR_B == "chr11") %>%
  mutate(
    BP_A_bin = cut(BP_A, breaks = seq(min(BP_A), max(BP_A), by = 50000), labels = FALSE, right = FALSE), # 50kb bins
    BP_B_bin = cut(BP_B, breaks = seq(min(BP_B), max(BP_B), by = 50000), labels = FALSE, right = FALSE)  # 50kb bins
  ) %>%
   group_by(BP_A_bin, BP_B_bin) %>% # Group by both bins
  summarise(avg_R2 = mean(R2, na.rm = TRUE), .groups = "drop") 

chr11_data_binned <- chr11_data_binned %>%%>% # Average R2 in each bin
  mutate(
    BP_A_mid = seq(min(BP_A), max(BP_A), by = 50000)[BP_A_bin],  # Midpoint of BP_A bin
    BP_B_mid = seq(min(BP_B), max(BP_B), by = 50000)[BP_B_bin]   # Midpoint of BP_B bin
  )


# 2. Create the heatmap with binned data
ggplot(chr11_data, aes(x = BP_A, y = R2)) +
  geom_point(color= "#e53207", alpha=0.4, size=0.5) +
  #scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Linkage Disequilibrium (R²) on Chromosome 11",
       x = "Base Position A",
       y = "R2",
       fill = "Average R²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################USE RGT_PST_FST_focal_regions to load in fst on chr 11

fst.chr11 <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/chr11.test.txt", header=F)
colnames(fst.chr11) <- c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2", "avg_wc_fst", "no_snps")

fst.chr11 <- fst.chr11 %>%
  filter(fst_pop_key == c("rustica.eu_rustica.r"))

ggplot(fst.chr11, aes(x = window_pos_1, y = avg_wc_fst)) +
  geom_point(color= "grey22", alpha=0.4, size=0.5) +
  #scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Linkage Disequilibrium (R²) on Chromosome 11",
       x = "Base Position A",
       y = "R2",
       fill = "Average R²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#https://ood-2.ai.lrz.de/rnode/cpu-008.ai.lrz.de/8933/graphics/plot_zoom_png?width=1920&height=527
####Look ino the SNP-Meth comparisons
#Took only the ones with a distance of 2 to begin with, can be done for all up to 500 or so base-pairs
snp.info <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/3_LD_decay/1_SNP_METH/4_HR_Rall_snp_meth_dist2.20.txt", header=T)
meth.info <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/prop.168.merged.20miss.rmlowv.txt", header=T)

#Remove depth information from SNPs
snp.info <- snp.info %>%
  mutate(across(contains("_ADL_"), ~ sub(":.*", "", .)))

#Create columns from which to join
snp.info <- snp.info %>%
  mutate(pos_minus = (POS - 2))%>%
  mutate(pos_plus = (POS + 2))

snp.info <- snp.info %>%
  mutate(site.id.m = paste0(CHROM,"_",pos_minus))%>%
  mutate(site.id.p = paste0(CHROM,"_",pos_plus))

joined.m <- inner_join(meth.info, snp.info, by=c("site.id" = "site.id.m"))
joined.p <- inner_join(meth.info, snp.info, by=c("site.id" = "site.id.p"))

joined <- rbind(joined.m[,1:350], joined.p[1:350])

joined_long <- joined %>%
  pivot_longer(
    cols = 2:169,  # Adjust the column range as needed
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  mutate(
    DataType = case_when(
      grepl("^ME_", Sample) ~ "METH_value",
      TRUE ~ "SNP_value"
    )
  ) %>%
  pivot_wider(names_from = DataType, values_from = Value)


joined_long <- joined_long %>%
  pivot_longer(
    cols = 13:180,  # Adjust the column range as needed
    names_to = "Sample_ID",
    values_to = "METH_value"
  )


joined_long <- joined_long %>%
  transform(Sample_ID=str_replace(Sample_ID,"ME_",""))

joined_long <- joined_long %>%
  filter(Sample_ID == Sample)

joined_long <- joined_long %>%
  transform(SNP_value=str_replace(SNP_value,"0/0","0|0"))
joined_long <- joined_long %>%
  transform(SNP_value=str_replace(SNP_value,"0/1","0|1"))
joined_long <- joined_long %>%
  transform(SNP_value=str_replace(SNP_value,"1/1","1|1"))

joined_long <- joined_long %>%
  filter(!SNP_value == "./.")

subset <- joined_long[c(1:12460),]

ggplot(subset, aes(x=as.factor(SNP_value), y=as.numeric(METH_value)))+
  geom_boxplot(color="blue") +
  facet_wrap(~ site.id) +
  theme_minimal()


###Check for farther distances
#Took only the ones with a distance of 2 to begin with, can be done for all up to 500 or so base-pairs
snp.info <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/3_LD_decay/1_SNP_METH/4_HR_Rall_snp_meth_dist2.20.txt", header=T)
snp.distance <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/2_Divergence/3_LD_decay/1_SNP_METH/3_HR_Rall_snp_meth_dist2.20.txt", header=T)
meth.info <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/168.20miss/prop.168.merged.20miss.rmlowv.txt", header=T)

#Remove depth information from SNPs
snp.info <- snp.info %>%
  mutate(across(contains("_ADL_"), ~ sub(":.*", "", .)))

#join to the distances table
snp.a <- snp.info %>%
      inner_join(snp.distance, by = c("POS" = "BP_A"), keep=T)
snp.b <- snp.info %>%
      inner_join(snp.distance, by = c("POS" = "BP_B"), keep=T)

snp.info <- bind_rows(snp.a, snp.b)

#Create columns from which to join
snp.info <- snp.info %>%
  mutate(pos_minus = (POS - distance))%>%
  mutate(pos_plus = (POS + distance))

snp.info <- snp.info %>%
  mutate(site.id.m = paste0(CHROM,"_",pos_minus))%>%
  mutate(site.id.p = paste0(CHROM,"_",pos_plus))

snp.info.long <- snp.info %>%
  pivot_longer(
    cols = 10:177,  # Adjust the column range as needed
    names_to = "Sample",
    values_to = "SNP"
  ) 

meth.info.long <- meth.info %>%
  pivot_longer(
    cols = 2:169,  # Adjust the column range as needed
    names_to = "Sample",
    values_to = "Meth"
  ) 

meth.info.long <- meth.info.long %>%
  transform(Sample=str_replace(Sample,"ME_",""))

joined.m <- inner_join(snp.info.long, meth.info.long, by=c("site.id.m" = "site.id", "Sample"), keep=T)
joined.p <- inner_join(snp.info.long, meth.info.long, by=c("site.id.p" = "site.id", "Sample"), keep=T)

joined <- rbind(joined.m, joined.p)

joined <- joined %>%
  transform(SNP=str_replace(SNP,"0/0","0|0"))
joined <- joined %>%
  transform(SNP=str_replace(SNP,"0/1","0|1"))
joined <- joined %>%
  transform(SNP=str_replace(SNP,"1/1","1|1"))

joined <- joined %>%
  filter(!SNP == "./.")

subset <- joined %>%
  filter(R2 > 0.7)

ggplot(subset, aes(x=as.factor(SNP), y=as.numeric(Meth)))+
  geom_boxplot(color="blue") +
  facet_wrap(~ site.id) +
  theme_minimal()

##other plots
pst_homer_ld <- inner_join(joined, pst_homer, by = c("site.id" = "site.id"))
pst_homer_ld_highr2 <- pst_homer_ld %>%
  filter(R2 > 0.5)
pst_homer_ld_highr2 <- pst_homer_ld_highr2 %>%
  group_by(site.id) %>%
  slice(1) %>%
  ungroup

pst_homer_ld <- pst_homer_ld %>%
  group_by(site.id) %>%
  slice(1) %>%
  ungroup

# Calculate category proportions
category_counts <- pst_homer_ld_highr2 %>%
  dplyr::count(category)

# Create pie chart
total <- ggplot(category_counts, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High R2 Categories") +
  theme_void()

####COMPARED TO WGS DISTRIBUTIONS
category_counts_homer <- pst_homer_ld %>%
  dplyr::count(category)

# Create pie chart
total_wgs <- ggplot(category_counts_homer, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Total Category") +
  theme_void()

total + total_wgs

###JUST region and annotation
pst_homer_ld <- pst_homer_ld %>%
  separate(Annotation, into = c("sea", "region"), sep = "_", remove = FALSE)
pst_homer <- pst_homer %>%
  separate(Annotation, into = c("sea", "region"), sep = "_", remove = FALSE)

# Count occurrences of `sea`
sea_counts <- pst_homer_ld %>%
  dplyr::count(sea)

# Count occurrences of `sea`
sea_wgs_counts <- pst_homer %>%
  dplyr::count(sea)

# Pie Chart for `sea`
p2 <- ggplot(sea_counts, aes(x = "", y = n, fill = sea)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "High LD CpG Regions") +
  theme_void()

p3 <- ggplot(sea_wgs_counts, aes(x = "", y = n, fill = sea)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Total CpG regions") +
  theme_void()

p2 + p3

# Count occurrences of `region`
region_counts <- pst_homer_ld %>%
  dplyr::count(region)
region_wgs_counts <- pst_homer %>%
  dplyr::count(region)

# Pie Chart for `region`
p4 <- ggplot(region_counts, aes(x = "", y = n, fill = region)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "High LD Region") +
  scale_fill_brewer(palette = "Dark2") + # Or scale_fill_manual(), scale_fill_viridis(), etc.
  theme_void() +
  theme(legend.position = "none")

p4 <- p4 + theme(legend_position ="none")

p5 <- ggplot(region_wgs_counts, aes(x = "", y = n, fill = region)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Total Regions") +
  theme_void()

p4+p5


##Find high LD region in chromsome11
chr11.ld <- snp.meth.ld %>%
  filter(CHR_A == 11)

chr11.ld <- chr11.ld %>%
  filter(R2 > 0.5)

##Find high LD region in chromsome11
chr03.ld <- snp.meth.ld %>%
  filter(CHR_A == 3)

chr03.ld <- chr03.ld %>%
  filter(R2 > 0.5)
