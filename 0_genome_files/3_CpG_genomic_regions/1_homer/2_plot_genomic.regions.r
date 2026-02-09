##Oct.2022
#Sarah Mueller
#script to evaluate missing data in methylation datasets
#Libraries
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
library(data.table)
library(methylKit)
library(gghalves)
library(wesanderson)

#sample names
sample_names <- read.table("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.168.txt")

#Dataset 1: 
#import homer output
homer.output.168.0miss <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/1_homer/homer.output.168.20miss.txt")
names(homer.output.168.0miss)[names(homer.output.168.0miss) == "PeakID..cmd.annotatePeaks.pl..dss.dsshome1.lxc0B.ra52qed.scripts.0.genome.files.homer.homer.input.168.10miss.bed..dss.dsslegfs01.pr53da.pr53da.dss.0034.assemblies.Hirundo.rustica.genome.v2.Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi.C.DLS_v2.fasta..gtf..dss.dsslegfs01.pr53da.pr53da.dss.0034.assemblies.Hirundo.rustica.genome.v2.annotation.GCF_015227805.1_bHirRus1.pri.v2_genomic.gtf."] <- "PeakID"
#separate annotation column to have separate bins
homer.output.168.0miss$Annotation = sub("\\(.*", "", homer.output.168.0miss$Annotation)
#convert . to underscore
homer.output.168.0miss$PeakID = sub("\\.", "_", homer.output.168.0miss$PeakID)
#convert annotation to factor
homer.output.168.0miss$Annotation <- as.factor(homer.output.168.0miss$Annotation)
summary(homer.output.168.0miss$Annotation)

#load in methylation data
load("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/2_methylkit.unite.filter/168.20miss/prop.168.20miss.RData")
joined.168.0miss <- left_join(prop.168.0miss, homer.output.168.0miss, by = c("site.id" = "PeakID"))
summary(joined.168.0miss$Annotation)
joined.168.0miss <- pivot_longer(joined.168.0miss, 2:169)
##plot proportion versus genomic region
n=0.05
pal <- c(wes_palette("Darjeeling1"))

ggplot(joined.168.0miss, aes(Annotation, value, fill=Annotation)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, 
                    width=0.5, nudge=n) +
  scale_fill_manual(values=pal)+
  geom_half_violin(side="r", nudge=n, alpha=0.3) +
  theme_minimal()

#Dataset2
#import homer output
homer.output.168.10miss <- read.delim("~/output/0.genome.files/3.homer/homer.output.168.10miss.txt")
names(homer.output.168.10miss)[names(homer.output.168.10miss) == "PeakID..cmd.annotatePeaks.pl..dss.dsshome1.lxc0B.ra52qed.scripts.0.genome.files.homer.homer.input.168.10miss.bed..dss.dsslegfs01.pr53da.pr53da.dss.0034.assemblies.Hirundo.rustica.genome.v2.Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi.C.DLS_v2.fasta..gtf..dss.dsslegfs01.pr53da.pr53da.dss.0034.assemblies.Hirundo.rustica.genome.v2.annotation.GCF_015227805.1_bHirRus1.pri.v2_genomic.gtf."] <- "PeakID"
#separate annotation column to have separate bins
homer.output.168.10miss$Annotation = sub("\\(.*", "", homer.output.168.10miss$Annotation)
#convert . to underscore
homer.output.168.10miss$PeakID = sub("\\.", "_", homer.output.168.10miss$PeakID)
#convert annotation to factor
homer.output.168.10miss$Annotation <- as.factor(homer.output.168.10miss$Annotation)
summary(homer.output.168.10miss$Annotation)

#load methylation data
load("./../4.0.methyl.unite.filter/prop.168.10miss.RData")
joined.168.10miss <- left_join(prop.168.10miss, homer.output.168.10miss, by = c("site.id" = "PeakID"))
summary(joined.168.10miss$Annotation)
joined.168.10miss <- pivot_longer(joined.168.10miss, 2:169)
##plot proportion versus genomic region
n=0.05
pal <- c(wes_palette("Darjeeling1"))

ggplot(joined.168.10miss, aes(Annotation, value, fill=Annotation)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, 
                    width=0.5, nudge=n) +
  scale_fill_manual(values=pal)+
  geom_half_violin(side="r", nudge=n, alpha=0.3) +
  theme_minimal()


###Dataset 3
#import homer output
homer.output.168.20miss <- read.delim("~/output/0.genome.files/3.homer/homer.output.168.20miss.txt")
names(homer.output.168.20miss)[names(homer.output.168.20miss) == "PeakID..cmd.annotatePeaks.pl..dss.dsshome1.lxc0B.ra52qed.scripts.0.genome.files.homer.homer.input.168.10miss.bed..dss.dsslegfs01.pr53da.pr53da.dss.0034.assemblies.Hirundo.rustica.genome.v2.Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi.C.DLS_v2.fasta..gtf..dss.dsslegfs01.pr53da.pr53da.dss.0034.assemblies.Hirundo.rustica.genome.v2.annotation.GCF_015227805.1_bHirRus1.pri.v2_genomic.gtf."] <- "PeakID"
#separate annotation column to have separate bins
homer.output.168.20miss$Annotation = sub("\\(.*", "", homer.output.168.20miss$Annotation)
#convert . to underscore
homer.output.168.20miss$PeakID = sub("\\.", "_", homer.output.168.20miss$PeakID)
#convert annotation to factor
homer.output.168.10miss$Annotation <- as.factor(homer.output.168.10miss$Annotation)
summary(homer.output.168.10miss$Annotation)

#load methylation data
load("./../4.0.methylunite/prop.168.10miss.RData")
joined.168.10miss <- left_join(prop.168.10miss, homer.output.168.10miss, by = c("site.id" = "PeakID"))
summary(joined.168.10miss$Annotation)
joined.168.10miss <- pivot_longer(joined.168.10miss, 2:169)
##plot proportion versus genomic region
n=0.05
pal <- c(wes_palette("Darjeeling1"))

ggplot(joined.168.10miss, aes(Annotation, value, fill=Annotation)) + 
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, 
                    width=0.5, nudge=n) +
  scale_fill_manual(values=pal)+
  geom_half_violin(side="r", nudge=n, alpha=0.3) +
  theme_minimal()


##############
#############
##############
##extra stuff for CpGislands and repeat

#import repeat and CpG island dataset
intersect.output.168.0miss <- read.delim("~/output/0.genome.files/3.CpG.genomic.regions/intersect.output.168.0miss.bed", header = FALSE)
colnames(intersect.output.168.0miss) <- c("Chr","Start","End","Strand","PeakID","remove", "remove.1", "remove.2","remove.3", "Annotation")

#import all regions bed intersect to see overlap
intersect.all.regions.168.0miss <- read.delim("~/output/0.genome.files/3.CpG.genomic.regions/intersect.all.regions.168.0miss.bed", header = FALSE)
colnames(intersect.all.regions.168.0miss) <- c("Chr","Start","End","Strand","PeakID","chr.region", "start.region", "end.region","blank", "Annotation", "strand.region")

#join so we can check homer vs. bedtools
#so essentially in the annotation gtf file regions are marked as several different categories (it can be in a transcript, as well as exon and CDS)
#in homer, there is a way to edit the parseGTF.pl which assigns priority to the different rankings, however, as far as i can tell there is no way to keep the several annotations so you can see all of them.
#this can be done by joining this to the genomic regions bed, however, this creates several rows for each entry, could use pivot wider to make it smaller
join.regions.homer <- left_join(homer.output.168.0miss, intersect.all.regions.168.0miss, by = "PeakID")

#select columns that match between both to create one merged file
homer.output.168.0miss <- homer.output.168.0miss %>%
  dplyr::select(PeakID, Chr, Start, End, Strand, Annotation)
intersect.output.168.0miss <- intersect.output.168.0miss %>%
  dplyr::select(PeakID, Chr, Start, End, Strand, Annotation)
intersect.remove <- inner_join(homer.output.168.0miss, intersect.output.168.0miss, by = "PeakID")
intersect.output.168.0miss <- intersect.output.168.0miss %>% 
  filter (! duplicated(PeakID))

genomic.regions <- rbind(homer.output.168.0miss)
genomic.regions$PeakID = sub("\\.", "_", genomic.regions$PeakID)
genomic.regions <- genomic.regions %>% distinct(PeakID, .keep_all= TRUE)