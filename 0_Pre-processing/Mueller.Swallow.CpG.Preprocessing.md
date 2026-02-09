## Metadata

Here I keep the main metadata file: hirundo_WGS_sampling_metadata_10.10.22_sm.xslx

Additionally, I have the plots of the sampling distribution. I use R packages that can work with shapefiles, so I downloaded shapefiles from the internet of the swallow distribution...which comes from Birdlife, and then simple country outlines and other climate information.

I thien plot these shapefiles along with my samples using: 1_Fig_1_Conceptual_Sample_Overview.Rmd

# Swallow DNA methylation

## 1.Data Import

See Data Import.md for more details. Loaded files into rawdata/

## 2. Preparing Genome Files

I will prepare some genome based files (CpG islands, cut-sites, etc.), to make downstream processing easier.

All files can be found within the project folder under 2021_SwallowRRBS/0_genome_files

#### 2.1 Cut-sites

The following is based on a script from Justin Merundon.

```bash
##make a new document with the following 2 lines to mimic the fasta file format for the cutsite
vim mspI
>mspI
CCGG

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta

seqkit locate ${genome} --bed -f mspI > mspI_5p_seqkit.bed
seqkit subseq ${genome} --bed mspI_5p_seqkit.bed > mspI_5p_seqkit_check.fa
```

Can also use oligoMatch:

```bash
#!/bin/bash

#SBATCH -A snic2018-3-655
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 6:00:00

oligoMatch mspI ${genome} mspI_5p
```

Then extend each cutsite by 150bp in either direction (since we cannot know strand-specificity with a CCGG cutsite)

```bash
bedtools slop -i mspI_5p_seqkit.bed -g chromosomes.genome -b 150 > mspI_5p_oligo_150bp.bed
```

Additionally, I want to make a new file that has cutsites only pesent on major chromsomes because there are many small contigs which are likely un-informative that I will remove:

```bash
grep 'chr' mspI_5p_seqkit.bed | wc -l
##1819728
##this makes a file with only those cutsites starting with 'chr'
head -1819728 mspI_5p_seqkit.bed >> mspI_5p_chr_seqkit.bed
##Need to extend that by 150bp (which represents the read length from sequencing in either direction, since it is paired end, this would cover any CpG sites that would be picked up by sequencing)
bedtools slop -i mspI_5p_seqkit.bed -g chromosomes.genome -b 150 > mspI_5p_oligo_150bp.bed
```

This can be used in filtering the CpG file later.

#### 2.2 Chromosome Lists

We also need a list with all chromosomes and unplaced scaffolds:

```bash
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta

samtools faidx ${genome}

mkdir genome.files

#The following files only have the major chromosomes
grep 'chr' ${genome}.fai > ./genome.files/chromosomes.list
#Convert 1-based index to 0-based .bed format 
awk '{OFS="\t"}{print $1, "0", $2-1}' ./genome.files/chromosomes.list > ./genome.files/chromosomes.bed
#Create a genome length file for bedtools as well, called a 'genome' file 
awk '{OFS="\t"}{print $1, $2-1}' ./genome.files/chromosomes.list > ./genome.files/chromosomes.genome


##I made a file that also has the scaffolds and unplaced, just in case this needs to be used as well
grep 'chr' ${genome}.fai > ./genome.files/chromosomes.scaffold.list
grep 'scaffold' ${genome}.fai >> ./genome.files/chromosomes.scaffold.list
grep 'unplaced' ${genome}.fai > ./genome.files/chromosomes.scaffold.list
```

#### 2.3 CpG Island Track (14.01.22)

We can use Takai and Jones’ algorithm for finding CpG islands in genomes implemented in [TaJoCGI]([GitHub - lucasnell/TaJoCGI: Find CpG islands in genome](https://github.com/lucasnell/TaJoCGI)) to get information on CpG islands within the swallow genome. The algorithm uses a sliding window approach (200 bp windows) to search for regions where the GC content is above 0.55 and GCobs/GCexp is greater than 0.65. These are then marked within the genome.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=6:00:00
#SBATCH -J TaJoCGI

python ~/scripts/TaJoCGI/TaJoCGI.py -o ~/genome.files/Hirundo_rustica_CGI.bed /dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta
```

This resulted in Hirundo_rustica_CGI.bed (29393 lines)

Then we can use [makeCGI]([Model Based CpG Islands](http://www.haowulab.org/software/makeCGI/index.html)) to create another CpG island track, using R (4.0.5):

NOTE: the fasta file needs to be in a folder labeled rawdata, with the suffix .fa not .fasta

```r
#setworking directory
  setwd("/dss/dsshome1/lxc0B/ra52qed/scripts/makeCGI")
##if using for the first time:
#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("rtracklayer")
#install.packages("makeCGI.tar.gz", repos = NULL, type = "source")

#use for iterations
library(makeCGI)
#Set up default parameters:
  .CGIoptions=CGIoptions()
#If the inputs are fa files, users need to create a directory called "rawdata" under current working directory and save all fa files there, then speicify the input type and species name. For example,
.CGIoptions$rawdat.type="txt"
.CGIoptions$species="Hirundo_rustica"
#Start running:
  makeCGI(.CGIoptions)
```

This resulted in CGI-Hirundo_rustica.txt (71788 lines)

Now we want to merge the results, which we can do with a simple bash script:

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=5:00
#SBATCH -J mergeCGI


taxa=Hirundo_rustica
#Output raw makecgi File
makecgi=CGI-Hirundo_rustica.txt
#Minimum cpg island size 
makecgiwin=250

#Filter for windows > window size 
awk -v w=${makecgiwin} '$4 > w' ${makecgi} > ${taxa}_makecgi_${makecgiwin}.tmp

#Output first three columns (chr, start, end) 
awk '{OFS="\t"}{print $1, $2, $3}' ${taxa}_makecgi_${makecgiwin}.tmp > ${taxa}_makecgi_${makecgiwin}.bed

#remove temporary file
rm *.tmp

#delete first line in makecgi .bed because it only contains identifiers
sed -i '1d' ${taxa}_makecgi_${makecgiwin}.bed

#Find the consistent islands between tajo's algorithm and makecgi
bedtools intersect -a ${taxa}_makecgi_${makecgiwin}.bed -b ${taxa}_CGI.bed > ${taxa}_CGI_merged.bed

#makeCGI island count
wc -l Hirundo_rustica_makecgi_250.bed
43288
#Tajo CpG island count
wc -l Hirundo_rustica_CGI.bed
29393
#Intersected results
wc -l Hirundo_rustica_CGI_merged.bed
30631 
```

This gives us 30,631 CpG islands. In the chicken there are [20,224](http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-galGal3.txt) identified. One thought on the discrepancy is the number of unplaced chromosomes so chr2_random etc. which is not present in the chicken. So it is possible that some would be collapsed if they were really placed on the genome

#### 2.4 Repeat Track (18.01.22)

We can identify repeat regions and separate methylation calls in these regions, since we have less confidence in these calls.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00
#SBATCH -J repeatmasker

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fna
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/

RepeatMasker ${genome} -species chicken -nolow -pa 20 -gff -dir ${genomedir}
```

This creates new files in the genome directory like fna.masked

We want to use these output files to create a bed file with all of the masked regions, to be able to remove these from the bismark extraction later on.

```bash
GFF=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta.out.gff
awk '{print $1, $4, $5, $7, $10}' $GFF > repeats.tmp1
grep -E -- 'chr' repeats.tmp1 > repeats_chrom.tmp
awk '{OFS="\t"}{print $1, $2, $3, $4, $5}' repeats_chrom.tmp | sed 's/"//g' > repeats.tmp2
cp repeats.tmp2 HR_repeats.bed
```

This will be used to filter out repeats later on.

#### 2.5 SNP Track

See WGS Data.md to see how I filtered the vcf file we recieved from the Safran Lab (Boulder, CO). After filtering then we want to get a SNP track that can be used for C-T filtering, etc.

In this case, I wanted to grab SNPs from both the vcf file from the Safran lab which contains more individuals and therefore may cover more SNPs that occur at low frequencies. Additionally, I want to pull SNPs from the merged vcf file that includes the European H.r.rustica individuals where there may be population specific SNPs which I do not want to miss.

Therefore, I run the same script twice on two vcf files. Importantly, here we filter out C-T SNPs as well as G-A SNPs. For the G-A SNPs, we must take the SNP position - 1, since we have already merged the CpG sites using the plus strand and therefore methylation will only show up in the plus strand, but a G-A SNP on the complementary CpG site on the minus strand would also cause mis-called methylation.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=60:00
#SBATCH -J vcf2bed

#output all positions from vcf file
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'REF="C" & ALT="T" | REF="G" & ALT="A"' /dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/Safran_lab/hirundo.rustica.300.id.allsites.renamed.miss02.maf05.ingroup.vcf.gz | awk '{
  if ($3 == "G" && $4 == "A") {
    print $1 "\t" ($2 - 1)
  } else if ($3 == "C" && $4 == "T") {
    print $1 "\t" $2
  }
}' > HR.300.m02.m05.d2_30.tmp

#take first 2 columns and convert to 0-based format
awk '{print $1 "\t" $2 "\t" $2}' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/2_CT_filter/HR.300.m02.m05.d2_30.tmp >> /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/2_CT_filter/HR.300.m02.m05.d2_30.bed
```

And again for the other vcf with concatenating both files:

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=60:00
#SBATCH -J vcf2bed

#Extract C-T and G-A and subtract one position from the G-A SNPs since we merge before filtering so these will show up as bp-1
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'REF="C" & ALT="T" | REF="G" & ALT="A"' /dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/merged/invariant/HR.168.filtered.vcf.gz | awk '{
  if ($3 == "G" && $4 == "A") {
    print $1 "\t" ($2 - 1)
  } else if ($3 == "C" && $4 == "T") {
    print $1 "\t" $2
  }
}' > HR.168.filtered.tmp

#Add the third column to the text file to create a bed file
awk '{print $1 "\t" $2 "\t" $2}' HR.168.filtered.tmp > HR.168.filtered.bed

#Concatenate the two files to create one bed file for extraction
cat HR.168.filtered.bed HR.300.m02.m05.d2_30.bed > all_ct_ga_snps.bed
```

## 3. Trim and Map

#### 3.1 Quality check of raw sequencing runs

Ran fastqc on all raw input files, output in folder /ra52qed/fastqc

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=4:00:00
#SBATCH -J fastqc.raw.data


while read f; 
do fastqc $f;
done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/NAMES.txt


multiqc .
```

multiqc reports in /Notes/documentation/multiqc

everything looks more or less good, one sample seems to be throwing way off in comparison to all the others

#### 3.2 Sample Names List

```bash
ls *__R1__* >> sample.names.txt
cat sample.names.txt | awk -F '__R1__' '{print $1}'


#this removes everyrthing after the acession IDS, so technical replicates are still there, but all the other information is removed...can be editted to parse in other ways 
```

#### 3.3 Trimming

So I want to use [TrimGalore]([TrimGalore/Trim_Galore_User_Guide.md at master · FelixKrueger/TrimGalore · GitHub](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md))! to remove low quality bases, adapter sequences as well as C's that were involved in end repair, which could lead to mis-called methylated sites.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00
#SBATCH -J trim.galore

#-q removes low quality bases (<20) first for RRBS samples, then moves onto adapter removal
#--fastqc runs fastqc on the resulting .fastq file
#--gzip compresses output files
#--retain_unpaired keeps unpaired reads that one read was removed due to being too short after adapter removal/trimming
#--rrbs and --paired for our specific data
trim_galore --rrbs --paired --quality 25 --fastqc --fastqc_args "--outdir /dss/dsshome1/lxc0B/ra52qed/fastqc/trim.galore.2501" --gzip --retain_unpaired -o /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowGenomics/trim_galore --cores 4 *__R1__* *__R2__*
```

"Read 1 sequences with adapter contamination will be trimmed a further 2 bp from their 3' end, and Read 2 sequences will be trimmed by 2 bp from their 5' end to remove potential methylation-biased bases from the end-repair reaction
All Read 2 sequences will be trimmed by 2 bp from their 5' end to avoid poor qualities or biases (e.g. M-bias for BS-Seq applications)"

Because the cut site is at CGG the first bases are always C/T (depending on if it is methylated or not) followed by GG. The first site is a true methylation site, which can be mapped to the genome. Because this is a directional library, the first C on the 5' end will be a true methylation...however on the 3' end a C is added during A tailing. This was corrected for using the --rrbs function in trim_galore which trims 2 extra bases after the A, which should catch the last 'fake' C. See this [tutorial](https://www.bioinformatics.babraham.ac.uk/projects/bismark/RRBS_Guide.pdf) and [forum](https://www.biostars.org/p/9475211/) for more information.

![Screenshot from 20220131 153855png](file:///home/mueller/Downloads/R1_trimgalore_sequence_content.png)

All of the R2 reads begin with an A in the first base. At first, this seems like something is wrong, however, TrimGalore should have already trimmed 2 bases after the tailing A, which it has (GA) if we compare to the raw reads. This 100% 'A' composition is because after MspI cutting, the fragment will end in `CCG`. Now the first `C` is in non-CG context and as such is almost universally unmethylated. The middle`C` is the one that gets filled in during the end-repair procedure, and it get's filled in with unmethylated cytosines. This, after the end repair, fill-in reaction and bisulfite treatment, the fragment will always end in `TTG`. Thus, R2 will always start with `CAA` see [here](https://www.researchgate.net/post/Why_do_some_G_bases_read_as_A_bases_after_bisulfite_conversion_PCR_ampifcation_and_sequencing) for more information.

![Screenshot from 20220214 141444png](file:///home/mueller/Downloads/R2_trimgalore_sequence_content.png)

#### 3.4 Mapping

First, we have to prep the genome. This does the c->t conversions and g->a needed for mapping.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=20:00
#SBATCH -J bismark.genome

bismark_genome_preparation /dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fna
```

Now, we map the reads to the reference genome.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=70
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=99:00:00
#SBATCH -J bismark.map

workdir=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowGenomics/trim.galore
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2
output=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowGenomics/bismark.map

while read f;
do bismark --parallel 5 --output_dir ${output} --genome ${genome} -1 ${f}__R1__*val_1* -2 ${f}__R2__*val_2*;
done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/sample.names.txt
```

Next, I need to sort the .bam files

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00
#SBATCH -J sort.bam

while read f;

        do samtools sort ${f}*.bam -o ${f}.sorted.bam;

done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/sample.names.txt
```

See [here]([Bismark PE/SE directional and non-directional alignments · Issue #208 · FelixKrueger/Bismark · GitHub](https://github.com/FelixKrueger/Bismark/issues/208)) if we get unusally low mapping rates for data as it is a very similar experiment, the suggestions are:

1. see if the different trimming mode has an influence on the stats (doublechecked trimming, should be fine)
  
2. relax the `--score_min` function somewhat for the paired-end alignments, e.g. to `--score_min L,0,-0.4`
  
3. you could also think of specifying `--unmapped` during the PE alignment steps, and then use the unmapped reads in single-end mode to get the most out of your data. The commands for the SE data would then be: R1,`bismark --genome /your/genome/ R1_unmapped.fq.gz` and R2, `bismark --genome /your/genome/ --pbat R2_unmapped.fq.gz`
  

#### 3.5 Mapping Efficiency Check

The following can quickly check the mapping efficiency of this step using the following script, which pulls from the .report file produced by bismark. Qualimap does not perform accurately on bismark results.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=5:00
#SBATCH -J bismark.map

while read f;

        do echo ${f}; awk 'NR==9 {print}' ${f}*report*;
        done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/sample.names.full.txt
~                                                                                       
```

## 4. Conversion Efficiency

Determining [conversion efficiency]([Controlling for bisulfite conversion efficiency with a 1% Lamda spike-in - Enseqlopedia](http://enseqlopedia.com/2016/10/controlling-for-bisulfite-conversion-efficiency-with-a-1-lamda-spike-in/)) using DNA that is known to be unmethylated (like mtDNA, although there has since been critic of this method or what is now more common using unmethylated lambda phage DNA) is an important step as it tells us the conversion rate of Cs to Ts in our samples. This can be done by running essentially the same pipeline that we will use on our samples, except applying this to the lamda phage reference genome.

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=70
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=36:00:00
#SBATCH -J bismark.map

workdir=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowGenomics/trim.galore
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Escherichia.phage.Lambda
output=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowGenomics/mapping.efficency


#bismark_genome_preparation ${genome}

while read f;
do bismark --parallel 5  --output_dir ${output} --genome ${genome} -1 ${f}__R1__*val_1* -2 ${f}__R2__*val_2*; 
echo ${f}; awk 'NR==9 {print}' ${f}*report*;
done < /dss/dsshome1/lxc0B/ra52qed/scripts/mueller/sample.names.full.txt

cd ${output}
bismark_methylation_extractor -p --no_overlap --bedgraph --cytosine_report --output ${output}/methylation.extraction --parallel 8 --buffer_size 20G --genome_folder ${genome} *bismark_bt2_pe.bam

multiqc .
```

After then running multiqc on the output of this script, the output below was generated.

![conversionefficiencymultiqcpng](file:///home/mueller/Downloads/conversion.efficiency.multiqc.png)

The percent of methylated cyotines (2nd column above), represents the inverse of the conversion efficieny, therefore our conversion efficiency is between 99.7-99.8%.

## 5. Methylation Extraction

I want to use the bismark_extraction to pull out the sites with CpG methylation and get the files that include the positions, and coverage of methylated/unmethylated reads.

#### 5.1 Bismark Extraction

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=99:00:00
#SBATCH -J bismark.extraction


#p is for paired end
##no_overlap is to only read CpG sites from overlapping R1 and R2 once

bismark_methylation_extractor -p --gzip --no_overlap --bedgraph --cytosine_report --output /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowGenomics/bismark.map/1_output_files --parallel 8 --buffer_size 20G --genome_folder /dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/ *bismark_bt2_pe.bam
```

#### 5.2 Bismark Extraction Output Files

The main file I will use from the output is the CpG_report for merging symmetric sites, however, let's review what the output from each of the Bismark_extraction function is and contains (all examples come from <u>CH_18_R_18CH024_BL_ADL_M__SRR17643018__R1__d26a809a141ff3606bbb012b518a4454)</u>

The input was the {Sample_Name}_val_1_bismark_bt2_pe.bam

1. **{Sample_Name}_val_1_bismark_bt2_pe.bedGraph.gz**
  
  Chromosome/START/END/Percent_Methylation
  

```bash
scaffold318    2701    2702    0
scaffold318    2704    2705    100
scaffold318    2713    2714    100
scaffold496    17553    17554    100
scaffold496    17670    17671    100
scaffold517    3945    3946    100
scaffold568    7463    7464    66.6666666666667
scaffold568    7466    7467    100
scaffold568    7470    7471    66.6666666666667
```

This is actually covered in the bismark.cov (see below), so this may not be a necessary file to keep.

2. **{Sample_Name}_val_1_bismark_bt2_pe.bismark.cov.gz**

 Chrom/Start/End/Percent_Methylation/Count_Methylated/Count_Unmethylated

```bash
scaffold318    2702    2702    0    0    5
scaffold318    2705    2705    100    5    0
scaffold318    2714    2714    100    5    0
scaffold496    17554    17554    100    1    0
scaffold496    17671    17671    100    1    0
scaffold517    3946    3946    100    1    0
scaffold568    7464    7464    66.6666666666667    2    1
scaffold568    7467    7467    100    3    0
scaffold568    7471    7471    66.6666666666667    2    1
scaffold568    7555    7555    100    3    0
```

3. **{Sample_Name}_val_1_bismark_bt2_pe.CpG_report.txt**

This is genome based, and shows every single 'C' positions in a CpG context, therefore there will be many sites that show 0/0 for methylated/unmethylated counts.

Chromsome/Position/Strand/Count_Methylated/Count_Unmethylated/Dinucleotide/Trinucleotide

```bash
scaffold318    5    +    0    0    CG    CGG
scaffold318    6    -    0    0    CG    CGT
scaffold318    26    +    0    0    CG    CGT
scaffold318    27    -    0    0    CG    CGT
scaffold318    410    +    0    0    CG    CGT
scaffold318    411    -    0    0    CG    CGT
scaffold318    416    +    0    0    CG    CGT
scaffold318    417    -    0    0    CG    CGC
scaffold318    424    +    0    0    CG    CGA
scaffold318    425    -    0    0    CG    CGT
```

4. **{Sample_Name}_val_1_bismark_bt2_pe.cytosine_context_summary.txt**

This is a summary of the methylation within different genomic contexts. A short document (65 lines for each sample)

```bash
upstream    C-context    full context    count methylated    count unmethylated    percent methylation
A    CAA    ACAA    176    220    44.44
C    CAA    CCAA    1917    1967    49.36
G    CAA    GCAA    52    425    10.90
T    CAA    TCAA    26    103    20.16
A    CAC    ACAC    520    801    39.36
C    CAC    CCAC    289    615    31.97
G    CAC    GCAC    173    416    29.37
T    CAC    TCAC    73    161    31.20
A    CAG    ACAG    1260    2653    32.20
C    CAG    CCAG    1684    9830    14.63
G    CAG    GCAG    3001    5217    36.52
T    CAG    TCAG    1585    3361    32.05
A    CAT    ACAT    952    921    50.83
C    CAT    CCAT    228    248    47.90
G    CAT    GCAT    40    204    16.39
T    CAT    TCAT    19    147    11.45
A    CCA    ACCA    39    81    32.50
C    CCA    CCCA    166    812    16.97
G    CCA    GCCA    78    84    48.15
T    CCA    TCCA    35    244    12.54
A    CCC    ACCC    50    186    21.19
C    CCC    CCCC    1216    3816    24.17
G    CCC    GCCC    163    699    18.91
T    CCC    TCCC    367    1506    19.59
A    CCG    ACCG    880    1082    44.85
C    CCG    CCCG    7363    30387    19.50
G    CCG    GCCG    1303    3215    28.84
T    CCG    TCCG    3274    8910    26.87
A    CCT    ACCT    18    31    36.73
C    CCT    CCCT    192    3324    5.46
G    CCT    GCCT    138    350    28.28
T    CCT    TCCT    136    2653    4.88
A    CGA    ACGA    1279459    1882498    40.46
C    CGA    CCGA    2719767    4658460    36.86
G    CGA    GCGA    2191149    3564257    38.07
T    CGA    TCGA    780837    1308668    37.37
A    CGC    ACGC    2023958    2910553    41.02
C    CGC    CCGC    4807857    11995973    28.61
G    CGC    GCGC    4563570    8610253    34.64
T    CGC    TCGC    1710654    3215577    34.73
A    CGG    ACGG    3400661    5025639    40.36
C    CGG    CCGG    8614706    10536213    44.98
G    CGG    GCGG    6237028    13256104    32.00
T    CGG    TCGG    3008301    4862249    38.22
A    CGT    ACGT    1290471    1808526    41.64
C    CGT    CCGT    2733563    4674601    36.90
G    CGT    GCGT    2190277    3041389    41.87
T    CGT    TCGT    1101837    1789570    38.11
A    CTA    ACTA    8    10    44.44
C    CTA    CCTA    2    57    3.39
G    CTA    GCTA    21    29    42.00
T    CTA    TCTA    4    14    22.22
A    CTC    ACTC    25    73    25.51
C    CTC    CCTC    255    1488    14.63
G    CTC    GCTC    91    349    20.68
T    CTC    TCTC    234    728    24.32
A    CTG    ACTG    552    1663    24.92
C    CTG    CCTG    2371    63783    3.58
G    CTG    GCTG    822    5099    13.88
T    CTG    TCTG    367    26976    1.34
A    CTT    ACTT    18    24    42.86
C    CTT    CCTT    32    1071    2.90
G    CTT    GCTT    42    150    21.88
T    CTT    TCTT    51    971    4.99
```

5. **{Sample_Name}_val_1_bismark_bt2_pe.M-bias.txt**
  
  File that shows across all reads from this individual, the percent methylation at each position of the read 1-150bp (because our reads are 150bp long maximum). This can show if there is some sort of technical bias at certain positions along the read...for example in RRBS sequencing you expect some non-random distribution at the beginning of reads due to cut-sites.
  
  ```bash
  CpG context (R1)
  ================
  position    count methylated    count unmethylated    % methylation    coverage
  1    8486092    8220754    50.79    16706846
  2    2713    3123    46.49    5836
  3    237    590    28.66    827
  4    336505    696438    32.58    1032943
  5    261273    510935    33.83    772208
  6    406506    609773    40.00    1016279
  7    403851    768120    34.46    1171971
  ```
  
6. **{Sample_Name}_val_1_bismark_bt2_PE_report.txt**
  

This is a summary of the mapping/alignment

```basic
Final Alignment report
======================
Sequence pairs analysed in total:    27881378
Number of paired-end alignments with a unique best hit:    18271283
Mapping efficiency:    65.5% 
Sequence pairs with no alignments under any condition:    6412839
Sequence pairs did not map uniquely:    3197256
Sequence pairs which were discarded because genomic sequence could not be extracted:    477

Number of sequence pairs with unique best (first) alignment came from the bowtie output:
CT/GA/CT:    9124899    ((converted) top strand)
GA/CT/CT:    0    (complementary to (converted) top strand)
GA/CT/GA:    0    (complementary to (converted) bottom strand)
CT/GA/GA:    9145907    ((converted) bottom strand)

Number of alignments to (merely theoretical) complementary strands being rejected in total:    0

Final Cytosine Methylation Report
=================================
Total number of C's analysed:    1279952049

Total methylated C's in CpG context:    90198288
Total methylated C's in CHG context:    2351163
Total methylated C's in CHH context:    3264316
Total methylated C's in Unknown context:    149578
```

7. **{Sample_Name}_val_1_bismark_bt2_pe_splitting_report.txt**
  
  A methylation report with some summary statistics on the called methylation
  
  ```basic
  Processed 18270806 lines in total
  Total number of methylation call strings processed: 36541612
  
  Final Cytosine Methylation Report
  =================================
  Total number of C's analysed:    687387707
  
  Total methylated C's in CpG context:    48686467
  Total methylated C's in CHG context:    1344289
  Total methylated C's in CHH context:    1788069
  
  Total C to T conversions in CpG context:    83327692
  Total C to T conversions in CHG context:    193520897
  Total C to T conversions in CHH context:    358720293
  
  C methylated in CpG context:    36.9%
  C methylated in CHG context:    0.7%
  C methylated in CHH context:    0.5%
  ```
  

### 5.3 Merge symmtric CpG sites

Next, I want to merge the symmetrical CpG sites from the different strands together to increase the coverage at each site. This is because in animal genomes the methylation from both strands generally reflects eachother and can therefore be combined. I use a script from Robin Cristofari ([penguin tools github]([penguin-tools/fuse_cpg_sites.py at master · rcristofari/penguin-tools · GitHub](https://github.com/rcristofari/penguin-tools/blob/master/fuse_cpg_sites.py))) to merge CpG sites. This works on the CpG_report files and outputs a merged.cov file which has the same structure as the bismark.cov.gz files and can be piped directly into the filtering steps outlined below. I ran this script on all samples.

```bash
for i in ./0_Pre-processing/1_bismark_map/1_output_files/*CpG_report.txt.gz; do python ./0_Pre-processing/1_bismark_map/0_merge_CpG.py --input $i; done 
```

Output from this script shows that indeed almost all sites are symmetrically methylated and therefore merging is apprpriate. One example of the merged.report.txt is:

```bash
# Count and mean coverage per CpG site, by strand:

    plus    minus    both
n_CpG    416511    413108    2031685
cov    11.5479    11.8878    55.9807

# Strand-specific coverage for CpGs sequenced on both strand:

plus-strand:    27.99 X
minus-strand:    27.99 X

Average variation between strands:    0.01

71.01% of all CpG sites are sequenced on both strands
```

#### 5.4 Methylation Extraction Filtering

I want to filter for several different parameters, that mostly deal with the quality of the called CpG sites from bismark. To make sure the called CpG are 'true' sites, we want to keep only those on major chromosomes, filter for depth at each CpG site (get rid of low depth and high depth (mapped to the wrong region)), as well as only keep sites that are within the predicted cut-sites as predicted by seqkit. Additionall, we want to filter out C-T, G-A SNPs and DMRs identified to be related to sex. Importantly, bedtools subtract needs the position to only be one site in the bed file otherwise it will not work properly.

```bash
#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=30:00
#SBATCH -J bismark.extraction

#Sarah Mueller --based on script from Justin Merundon
#For filtering merged cov files from methylation extraction in bismark
#Usage: for i in $(ls *_val_1_bismark_bt2_pe.CpG_merged.cov | rev | cut -c37- | rev); do sbatch -J ${i} /dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/1_RRBS/0_Pre-processing/1_trim_map_filter/5_filter/1_bismark.filter.merged.sh ${i}; done


SCRATCH=/tmp
RUN=$SLURM_JOB_NAME


#Feature files 
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2
cgi=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/HR_CGI_merged.bed
snps=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/2_CT_filter/all_ct_ga_snps.bed
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/2_repeatmasker/HR.repeats.bed
chrs=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome
chrsbed=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.bed
sexdmr=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/3_sex_evaluation/sex_dmr_remove.bed
cutsites=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_cut_sites/mspI_5p_chr_seqkit_150bp.bed

#Directory paths 
workdir=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/1_bismark_map/covfiles
output=/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_output_files

#Need to copy the input files into the 'SCRATCH' directory (this could be the fastq file, bismark extraction, bam etc.) 
#Since this is for filtering, then its the bismark extraction output that we copy over 
cat ${workdir}/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov > $SCRATCH/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov

cd ${SCRATCH}

##Variables
#Minimum Coverage Required for a site to be kept
mincov=5
maxcov=200


#Only keep sites on the major chromosomes 
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' > ${SCRATCH}/${RUN}.tmp1

#Only keep positions within predicted inserts
bedtools intersect -wb -a ${cutsites} -b ${SCRATCH}/${RUN}.tmp1 | awk '{OFS="\t"}{print $7, $8, $9, $10, $11, $12}' | sort | uniq > ${SCRATCH}/${RUN}.tmp2

#Remove C-T SNP positions on major chromosomes
bedtools subtract -a ${SCRATCH}/${RUN}.tmp2 -b ${snps} > ${SCRATCH}/${RUN}.tmp3

#Remove sex related positions on major chromosomes
bedtools subtract -a ${SCRATCH}/${RUN}.tmp3 -b ${sexdmr} > ${SCRATCH}/${RUN}.tmp4

#Filter for sites with at least 5 reads and no more than 200 reads
awk -v x=${mincov} -v y=${maxcov} '($5 + $6) >= x && ($5 + $6) <= y' ${SCRATCH}/${RUN}.tmp4 > ${SCRATCH}/${RUN}.tmp5

#Filter out those from repeat regions 5mCs
bedtools intersect -wb -a ${repeats} -b ${SCRATCH}/${RUN}.tmp5 | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | gzip -c > ${output}/${RUN}.CpG_in_Repeats.cov.gz
bedtools subtract -a ${SCRATCH}/${RUN}.tmp5 -b ${repeats} | gzip -c > ${output}/${RUN}.CpG_5mC.cov.gz

### SUMMARIZE COUNTS
#Starting, raw positions
echo "RAW" >> ${output}/${RUN}.COUNTS.txt
cat ${workdir}/${RUN}_val_1_bismark_bt2_pe.CpG_merged.cov| wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining on major scaffolds
echo "MAJOR_CHROM" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp1 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after cut site filter
echo "CUT_SITE_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp2 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering SNPs
echo "SNP_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp3 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering sex related positions
echo "SEX_DMR_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp4 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions remaining after filtering for coverage
echo "COVERAGE_FILTER" >> ${output}/${RUN}.COUNTS.txt
cat ${SCRATCH}/${RUN}.tmp5 | wc -l >> ${output}/${RUN}.COUNTS.txt

#Positions in repeats
echo "REPEAT_FILTER" >> ${output}/${RUN}.COUNTS.txt
zcat ${output}/${RUN}.CpG_5mC.cov.gz | wc -l >> ${output}/${RUN}.COUNTS.txt
```

Run in each sample:

Note: the script must be slightly updated based on different naming structures (if youre file ends differently it will need to be updated), so the .cov.gz needs to be updated if this is not the proper ending.

```bash
for i in $(ls *_val_1_bismark_bt2_pe.CpG_merged.cov | rev | cut -c37- | rev); do sbatch -J ${i} /dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/1_RRBS/0_Pre-processing/1_trim_map_filter/5_filter/bismark.filter.sh ${i}; done
```

#### 5.5 Quality Checking

I use a script that is located within

```bash
./0_Pre-processing/2_bismark_filter/4_bismark_filter_plot.r
```

in order to look at the quality of mapping and extraction. This pulls from several files that are not within the bismark_filter folder but from multiqc, etc. so some paths may not work.

First, I need to compile the count files into a workable format for R, ideally with the sample name included and the data horizontally instead of vertically oriented:

```bash
for i in *COUNTS*; do tr "\n" "\t" < $i > $i".sorted"; done
rm *COUNTS.txt
awk '{print FILENAME" \"" $0"\""; nextfile}' *COUNTS* >> a.COUNTS.merged.tmp
sed -i 's/"//g' a.COUNTS.merged.tmp
awk '{print $1, $3, $5, $7, $9, $11}' a.COUNTS.merged.tmp >> a.COUNTS.clean.tmp
sed -i '1 i\Sample_ID\tRAW\tMAJOR_CHROM\tCUT_SITE_FILTER\tSNP_FILTER\tCOVERAGE_FILTER' a.COUNTS.cleaned.txt
awk -F '__' '$1=$1' OFS="\t" a.COUNTS.clean.tmp >> a.COUNTS.clean.tmp1
awk '{ print $1, $5, $6, $7, $8, $9}' a.COUNTS.clean.tmp1 >> a.COUNTS.clean.txt

rm *tmp*
```

##### 5.5.1 Percent of Aligned Reads per Sample

I plotted the Percentage of Aligned Reads per sample to see how well mapping and alignment performed. This can be especially important as RRBS with default settings can have low mapping rates, making parameter testing key (see above). This uses a 'cleaned general stats file' with shorter sample IDs (adapt script in next section to suit the needs of the sample names.)

```r
##map percent aligned from multiqc
df.multiqc <- read.table(file = "0_Pre-processing/0_fastqc/2_bismark.map.extraction/multiqc_general_stats.txt", header=T)
colnames(df.multiqc)<- c("Sample_ID","p_cpg","p_chg", "p_chh", "total_reads", "aligned_reads", "percent_aligned")

df.multiqc <- df.multiqc %>%
  separate(Sample_ID, into= c("Sample_ID", "SRR", "read", "id"), sep =  "__")

ggplot(df.multiqc,aes(x=Sample_ID,y=percent_aligned))+
  geom_bar(stat="identity",position = "dodge")+
  theme_classic(base_size=16)+
  ylim(0,100)+
  ggtitle("% Aligned Reads using Bismark")+
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())
```

##### 5.5.2 Filtering Output

Next, we want to plot the different filtering steps and corresponding counts into a histogram plot to make it easy to see how many CpG sites have been filtered out at each step:

```r
df <- read.delim(file = "./0_Pre-processing/2_bismark_filter/1_output_files/a.COUNTS.clean.txt", sep= " ", header=T)

#need to get data into a workable format
df.melt <- melt(df)

df.summary <- df.melt %>% 
  group_by(variable) %>%
  summarize(mean=mean(value),sd=sd(value))

#make color palette
pal <- c( wes_palette("Darjeeling1")[1],
          wes_palette("Cavalcanti1")[2],
          wes_palette("Darjeeling1")[5],
          wes_palette("Darjeeling2")[2],
          wes_palette("FantasticFox1")[2],
          wes_palette("FantasticFox1")[5],
          "grey33")

##make histogram
ggplot(df.melt,aes(x=SAMPLE_ID, fill=variable, y=value))+
  geom_bar(stat="identity",position="dodge")+
  theme_classic(base_size=18)+
  scale_fill_manual(values=pal)+
  ylab("Number of CpGs")+
  theme(axis.text.x=(element_text(angle=90,hjust=1)), axis.title.x=element_blank())


#Specifically, how many SNPs were filtered out
df <- df %>%
  mutate(snps_filtered = CUT_SITE_FILTER - SNP_FILTER)

summary(df$snps_filtered)
```

##### 5.5.3 M-Bias Plot

```bash
#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00
#SBATCH -J bismark.mbias

for RUN in $(ls *M-bias.txt | sed 's/\..*//g' ); do 

sed -n '/CHG context (R1)/q;p' ${RUN}.*M-bias.txt | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R1","CpG"}' > ${RUN}.MBIAS1.txt

sed -n '/CpG context (R2)/,$p' ${RUN}.*M-bias.txt | sed -n '/CHG context (R2)/q;p'  | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R2","CpG"}' > ${RUN}.MBIAS2.txt

done
echo -e 'Sample\tPosition\tCountMethylated\tCountUnmethylated\tPercentMethylated\tCoverage\tRead\tContext' > MBias.header
cat MBias.header *MBIAS1.txt > Methylation_Bias-Input_R1.txt
cat MBias.header *MBIAS2.txt > Methylation_Bias-Input_R2.txt

#clean-up
rm MBias.header
rm *M-Bias.txt*
```

And plot in R:

```r
#mbias plotting
bias <- read.table('Methylation_Bias-Input.txt',header=T)

a <- ggplot(bias,aes(x=Position,col=Sample))+
  geom_line(aes(y=PercentMethylated),stat="identity",show.legend=F,size=0.5)+ylab('Percent Methylation')+
  theme_classic(base_size=16)+
  scale_color_viridis(discrete=T,option='turbo')+
  coord_cartesian(ylim=c(0,100))+
  ggtitle('Methylation Along Read Length')

a

b <- ggplot(bias,aes(x=Position,col=Sample,y=log(Coverage)))+
  geom_line(stat="identity",show.legend=F,size=0.5)+ylab('log(Coverage)')+
  theme_classic(base_size=16)+
  scale_color_viridis(discrete=T,option='turbo')+
  ggtitle('Coverage Along Read Length')
b
```

At first, I thought the M-bias plots look weird (see below). But upon closer inspection, they closely reflect the multiqc trimming, which makes sense.

Forward:

![MbiasR1multiqcpng](file:///home/mueller/Downloads/M-bias_R1_multiqc.png)

For the forward (R1) reads the spike at the beginning reflects that the start position is always a CpG site, which means there is bias for the reasd position, but would not affect the overall results. Then there are smaller spikes along the read length, which is common for RRBS data and could result from i) repetitive regions, ii) mis-mapping, iii) conversion failure for some technical reason at specific CpG loci.

Reverse:

![MbiasR2multiqcpng](file:///home/mueller/Downloads/M-bias_R2_multiqc.png)

For the reverse reads, the first position having 0% methylation makes sense because all R2 at the 5' end will start with A, and since A cannot be methylated this will bias the reads (see trimming above).

The wider variation towards the 3' end and the unmethylated or methylated sites at this end usually are caused by the strict adapter trimming where even a single A is removed from the end, this will bias the last position, however, since not all reads are exactly 150bp (some are 147, 148, 149) this looks messy in graph form. For more information on 'weird' M-bias plots, [this]((http://seqanswers.com/forums/showthread.php?t=90956)) is a nice forum. If we take a subset, we can see the individual spikes for each sample:

![MethylationBiasR2subset1png](file:///home/mueller/Downloads/Methylation_Bias_R2_subset1.png)

#### 5.6 Coverage along Chromosomes (25.04.22)

Next, I want to explore the coverage along the chromsomes of called methlyation sites. I will go first with 1 sample to see if I can get a nice plot.

```bash
#total number of methylated cyotsines
for i in $(ls *.CpG_5mC.cov.gz | rev | cut -c16- | rev);do  zcat ${i}.CpG_5mC.cov.gz | awk -v s=${i} 'BEGIN{FS="\t"; sum=0} {sum+=$5} END{print s, sum}' >> coverage.txt ; done

#total number of unmethylated cytosines
for i in $(ls *.CpG_5mC.cov.gz | rev | cut -c16- | rev);do  zcat ${i}.CpG_5mC.cov.gz | awk -v s=${i} 'BEGIN{FS="\t"; sum=0} {sum+=$6} END{print s, sum}' >> coverage.txt ; done

#total number of cytosines (aka total number of rows)
for i in $(ls *.CpG_5mC.cov.gz | rev | cut -c16- | rev);do  zcat ${i}.CpG_5mC.cov.gz | awk -v s=${i} 'END{print s, NR}' >> total.lines.txt ; done

##combine these into one file with the columns Sample.ID, #methylated, # unmethylated, #total cytosine sites
##then we can plot them in R
```

Next, plotting in R:

```r
bs.genome <- toGRanges(data.frame(chr=c("chr01"), start=c(1), end=c(119023420)))

cov <- read.table("chr01.txt",header=FALSE)
head(cov)

cov_GR <- makeGRangesFromDataFrame(df=cov,ignore.strand=TRUE,seqinfo = NULL, seqnames.field = "V2",start.field = "V3", end.field = "V3",starts.in.df.are.0based=FALSE, keep.extra.columns=TRUE)
cov_GR

samples <- as.vector(unique(cov_GR$V1))
samples

total.tracks <- 28

#Plot
pp <- getDefaultPlotParams(plot.type=4)
pp$leftmargin <- 0.09
kp <- plotKaryotype(plot.type=4, genome = bs.genome, main = "chr01 coverage",plot.params = pp,labels.plotter=NULL)
kpAddChromosomeNames(kp,srt=90, cex=1)


for (i in 1:10) {
  name <- samples[i]
  sub <- cov_GR[grepl(name,cov_GR$V1), ]
  at <- autotrack(current.track = i, total.tracks = 10)
  kpPoints(kp,data=sub,y=sub$V5,cex=0.75,col="black",pch=16,r0=at$r0,r1=at$r1,ymin=0,ymax=200)
  kpAxis(kp,ymin=0,ymax=200,r0=at$r0,r1=at$r1)
}

#histogram of coverage between two chromsomes, or more depending on how you plot
#mean and std of chromosomes

chr01 <- cov %>% 
  group_by(V1,V2) %>%
  summarise(mean = mean(V5),n=n())


ggplot(chr01,aes(x=factor(V1),y=mean,fill=V2))+
  geom_bar(stat="identity",position = "dodge")+
  theme_classic(base_size=16)+
  ggtitle("Average Coverage across Chromosome 01 per Individual")+
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank(), legend.position="none")


coverage <- read.table("coverage.txt", header=T)

ggplot(coverage,aes(x=factor(Sample.Id),y=coverage,fill=wes_palette("Darjeeling1")[1]))+
  geom_bar(stat="identity",position = "dodge")+
  theme_classic(base_size=16)+
  ggtitle("Average coverage across all chromsomes per technical replicate")+
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank(), legend.position="none")
```

#### 5.7 Sex Identification/Confirmation against WGS (25.04.22)

I thought it would be a good idea to double check the sex of the individuals as this could point out abnormailities in sequencing, and a way to check against something we already know. To do this, we can pull the W and Z chromosomes out from the bismark.cov.gz file.

```bash
for i in $(ls *.CpG_5mC.cov.gz | rev | cut -c16- | rev);do 
zcat ${i}.CpG_5mC.cov.gz | egrep -w 'chrZ|chrW' | awk -v s=${i} '{OFS="\t"}{print s, $1, $2, $4, $5+$6}' >> Sex_ID.txt
done 
```

All:

![sexID1Zpng](file:///home/mueller/Downloads/sex.ID.1Z.png)

Males:

![sexidmpng](file:///home/mueller/Downloads/sex.id.m.png)

Females:

![sexidfpng](file:///home/mueller/Downloads/sex.id.f.png)

Unknown

![sexidupng](file:///home/mueller/Downloads/sex.id.u.png)

## 6. Evaluating Samples run on 2 lanes

Now, I want to evaluate how similar the samples that were run on two lanes are. These are samples from the same library but the first 'replicate' did not have enough coverage, so it was sequenced again. I want to double check that the different runs are highly correlated and there are no mistakes.

In order to do this, we can use methylkit to look at the correlation between technical replicates.

We will need a list of sample names and the sample with only SRIDs for R. We can add quotation marks to all sample files for R input:

```bash
awk '{ print "\""$0"\""","}' sample.list.txt
```

Next we can run in R:

(this is a subset of only a few samples, because the actual script is very long)

```r
##MethylKit Code for Swallow Epigenetics 
##began 02.03.2022
##Sarah Mueller

#load librraies
library(methylKit)

#set working directory
setwd("/home/mueller/Desktop/Notes/documentation/3.2.methylKit/")
list.files()

##lines from Justins script
##### METHYLKIT: List CpG cov reports #####
file.list=list("CH_18_R_18CH024_BL_ADL_M__SRR17643018__R1__d26a809a141ff3606bbb012b518a4454.CpG_5mC.cov.gz",
               "CH_18_R_18CH081_BL_ADL_M__SRR17643017__R1__4a1dc631af59269ba7ecf92fe45ea120.CpG_5mC.cov.gz")

CpG.bismark.covM=methRead(file.list, sample.id=list("CH_18_R_18CH024_BL_ADL_M__SRR17643018",
                                                    "CH_18_R_18CH081_BL_ADL_M__SRR17643017"),assembly="hirundo",treatment=c(1, 1), context ="CpG", pipeline= "bismarkCoverage")

CpG.bismark.covM

#Filter by coverage
filtered.CpG.bismark.covM=filterByCoverage(CpG.bismark.covM,lo.count=5,lo.perc=NULL,
                                           hi.count=NULL,hi.perc=95)

#Divide into technical replicates
tech1 <- reorganize(filtered.CpG.bismark.covM,sample.ids=c("CN_14_G_1089_BL_ADL_M__SRR17643192",
                                                           "CN_14_G_1089_BL_ADL_M__SRR17643274"),treatment=c(1,1))
meth_tech1=unite(tech1,mc.cores=4)
getCorrelation(meth_tech1,plot=FALSE,nrow=1000000)
t1 <- capture.output(getCorrelation(meth_tech1,plot=FALSE,nrow=1000000))
write.table(t1[3], file="t1.txt",col.names=FALSE,row.names=FALSE, append=T)

getCorrelation(meth_tech1,plot=TRUE)
png("t1.png",width=1200,height=1000)
dev.off()


tech2 <- reorganize(filtered.CpG.bismark.covM,sample.ids=c("CN_14_G_1097_BL_ADL_F__SRR17643047",
                                                           "CN_14_G_1097_BL_ADL_F__SRR17643058"),treatment=c(1,1))
meth_tech2=unite(tech2,mc.cores=4)
getCorrelation(meth_tech2,plot=FALSE,nrow=1000000)
t2 <- capture.output(getCorrelation(meth_tech2,plot=FALSE,nrow=1000000))
write.table(t2[3], file="t2.txt",col.names=FALSE,row.names=FALSE, append=T)
```

Ran this for all 86 replicates. Correlations were very high in all. We can visualize them in R:

```r
##plot technical replicates 
tech.reps <- read.table("reps.tech.txt", header=F)

pal <- c( wes_palette("Darjeeling1")[2],
          wes_palette("Cavalcanti1")[5])

ggplot(tech.reps, aes(x=factor(V1),y=V2))+
  geom_point(stat="identity", shape=21, size = 4)+
  scale_color_manual(values=pal)+
  theme_classic(base_size=16)+
  ggtitle("Technical Replicate Correlation")+
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank(), legend.position="none")
```

![6techrepscorrelationpng](file:///home/mueller/Downloads/6.tech.reps.correlation.png)

Then we want to merge all the bam files from the 2 or more 'technical replicates' into one 'merged bam file' from which we will make methylation calls again.

```bash
#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=10:00:00
#SBATCH -J bismark.mbias

for i in $(ls *.sorted.bam | rev | cut -c62- | rev); do samtools merge ${i}.merged.bam ${i}*.sorted.bam; done
```

Then run bismark.extraction again.

## 7. MethylKit Raw Data

Methylkit is an R package that can be used for DNA methylation data. It can do filtering, mostly based on coverage, normalization, and uniting. Uniting will merge the methylation data at sites along the genome from each individual. But first we need to look at the raw data and see what is going on there.

Methylkit can also perform count normalization, but only to a certain extent. If there are expected large batch effects, a better approach may be [ComBat]([Batch Correction | Griffith Lab](https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/)).

First, we can look at the distribution of methylation per sample, below is one example. If one sample looks very different than the rest, this may indicate issues with sequencing.

![7MethylkitrawdataCH008png](file:///home/mueller/Downloads/7.Methylkit.raw.data.CH008.png)

This means I have 600,000 sites of over 1 million which are unmethylated. Since we are using [RRBS we are capturing CpG islands which are usually unmethylated](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7655987/), so this type of distribution would make sense. Additionally, there are likely many sites which are non-variable across all individuals, and later we can remove these sites which are uninformative.

Coverage looks pretty decent, I would say over 20x is a very good start for methylation data:

![histogramcoverage946998png](file:///home/mueller/Downloads/histogram.coverage.946998.png)

Coverage across chromosome 1 in all individuals:

![chr01coveragepng](file:///home/mueller/Downloads/chr01.coverage.png)

There are clearly regions that have very high coverage (>=200). This could result from mis-mapping of repetitive regions. This can also be filtered out in methylkit by filtering for maximum coverage as well.

### 7.1 Creation of Dataset

For 85 samples from 6 different subspecies, I created one dataset to be used for downstream analysis. This dataset includes all 85 individuals, removes any sites above the 99% quantile, and unites methylation sites allowing up to 20% missing data per loci. For some downstream analyses, sites with missing data had to be excluded, this is documented in the methods.

The workflow is:

```r
library(methylKit)
library(readxl)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_output_covfiles")
list.files()

#CHANGE FOR EACH NEW METHYLOBJECT
#name of how you want the data to be saved
NAME= "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/meth.85.merged.20miss.txt"
N = 85
MIN.COV = 20
MAX.COV= 99.9


# Metadata for all samples
info <- read_xlsx("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/NCBI_info_updated_20220727.xlsx")

directory="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_output_covfiles"

filenames=list.files(path=directory, full.names=TRUE)
filenames=as.list(filenames)
names=list.files(path=directory)
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
meth=unite(normalized.myobj, min.per.group = 68L, destrand = F)

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
```

This unites all samples into one dataframe and removes misc. chromosomes/scaffolds. This outputs a matrix where there are 3 columns for each individual (methylated, unmethylated and total counts). In several analyses we cannot use the count of methylated Cs, but rather need proportion of methylation.

The following script (./0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/2_count.168.20miss.Rmd) creates those dataframes that will be used for all downstream analysis. This script also removes CpG loci that are hypo- or hyper- methylated (below 10% and above 90% in all individuals), removes low variability sites (loci that show the 5% smallest standard deviation across all individuals), and converts outlier methylation sites to NA.

Outlier methylation is an issue when most individuals cluster around 0.2 methylation, and 1-2 individuals have 0.8 methylation. These two individuals could be technical errors or true biological differences, however, either way they cause issues in downstream analysis. Therefore, here we have converted these values to NA. I also performed several analyses excluding these sites as well.

From there, I extract it to a dataframe and then create a counts file (this includes only the methylated read count) and then a table of proportion of methylated reads.

The sample names file just includes one sample name per line in the order they show up in in the methylkit object. I then add ME (stands for methylated), UN (stands for unmethylated), and TO (stands for total read count)

```r
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss")
#sample names
sample_names <- read.table("/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.85.txt")
#assign meth, total and unmethylated to sample names
sample_names <- cbind(sample_names, sample_names, sample_names)
colnames(sample_names) <- c("ME", "TO", "UN")
sample_names$ME <- sub("^", "ME_", sample_names$ME)
sample_names$UN <- sub("^", "UN_", sample_names$UN)
sample_names$TO <- sub("^", "TO_", sample_names$TO)

#sample size
N =85

##load methylation data
meth.85.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/meth.85.merged.20miss.txt", header=T)
```

I want to have a dataframe with the sample names and methylation data, as a proportion and counts of methylation.

```r
#grab methylation data
meth_basic_stat <- meth.85.20miss
#meth_basic_stat <- methylKit::getData(meth.168.0miss)
##get proportion data by taking number of methylated/total coverage
proportion_basic_stat <- (meth_basic_stat[ , c(unique(c(seq(6, ncol(meth_basic_stat), 3))))]/meth_basic_stat[ , c(unique(c(seq(5, ncol(meth_basic_stat), 3))))])
#change names to sample names for easy identification
colnames(proportion_basic_stat) <- sample_names[,1]
meth_basic_stat<- meth_basic_stat %>%
  tidyr::unite(site.id, chr:start)
rownames(proportion_basic_stat) <- meth_basic_stat[,1]

#move site.id to first column
proportion_basic_stat <- cbind(meth_basic_stat[,1], proportion_basic_stat)
colnames(proportion_basic_stat) <- c("site.id", sample_names[,1])
```

During this I will also remove sites that have a mean methylation of 10% or less and above 90%. After, I will calculate the standard deviation for the remaining sites, calculate a 5% threshold and remove all sites that show low variablility (standard deviation falls within the 5% threshold).Different papers use different thresholds for removing low variability sites (some 5%, 95%), so this is quite a loose threshold and could be further filtered later on, however, the order of filtering does matter. If you first filter for low variability,you are primarily capturing only hypo-methylated sites as most CpG islands are unmethylated.

```r
#calculate mean and standard deviation across individuals for all sites
proportion_basic_stat <- proportion_basic_stat %>%
  mutate(proportion.ave = rowMeans(proportion_basic_stat[,2:(N+1)], na.rm = TRUE),
    proportion.sd = apply(proportion_basic_stat[,2:(N+1)], 1, function(x) sd(x, na.rm = TRUE)))
#add site.id
#proportion_basic_stat<- proportion_basic_stat %>%
#                        mutate(site.id = meth_basic_stat[,1])

##total sites before filtering
ggplot(proportion_basic_stat, aes(x = proportion.ave)) +
  geom_histogram(binwidth = 0.05, fill = "#009193", color = "black", alpha = 0.7) +
  labs(title = "Average Methylation Proportion across Individuals Before Filtering",
       x = "Average Methylation Proportion",
       y = "Frequency") +
  theme_minimal(base_size = 12) +  # Set a base size for text
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
##before sd
ggplot(proportion_basic_stat, aes(x = proportion.sd)) +
  geom_histogram(binwidth = 0.05, fill = "#009193", color = "black", alpha = 0.7) +
  labs(title = "Standard Deviation across Individuals Before Filtering",
       x = "Standard Deviation",
       y = "Frequency") +
  theme_minimal(base_size = 12) +  # Set a base size for text
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
```

Now, filter out the hypo- and hyper- methylated sites

```r
#make df with filtered out sites, to run some plots on later
filtered_out <- proportion_basic_stat %>%                       
                        rowwise %>%
                        filter(proportion.ave < 0.1) %>%
                        filter(proportion.ave > 0.9)

#filter out hypo- (<y10%, hyper >90%,)
proportion_basic_stat<- proportion_basic_stat %>%                       
                        rowwise %>%
                        filter(proportion.ave > 0.1) %>%
                        filter(proportion.ave < 0.9)

#total sites after hypo- hyper- filtering
ggplot(proportion_basic_stat, aes(x = proportion.ave)) +
  geom_histogram(binwidth = 0.05, fill = "#009193", color = "black", alpha = 0.7) +
  labs(title = "Average Methylation Proportion across Individuals After Hypo/Hyper- Filter",
       x = "Average Methylation Proportion",
       y = "Frequency") +
  theme_minimal(base_size = 12) +  # Set a base size for text
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
##before sd
ggplot(proportion_basic_stat, aes(x = proportion.sd)) +
  geom_histogram(binwidth = 0.05, fill = "#009193", color = "black", alpha = 0.7) +
  labs(title = "Standard Deviation across Individuals After Hypo/Hyper- Filter",
       x = "Standard Deviation",
       y = "Frequency") +
  theme_minimal(base_size = 12) +  # Set a base size for text
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
```

```r
# Determine the 95th percentile threshold for the standard deviations
threshold <- quantile(proportion_basic_stat$proportion.sd, 0.05, na.rm = TRUE)

#add filtered out sites for plotting
filtered_out_lowv <- proportion_basic_stat %>%                       
                        rowwise %>%
                        filter(proportion.sd <= threshold)

# filter out low variability 5% threshold for standard deviation
proportion_basic_stat<- proportion_basic_stat %>%                       
                        rowwise %>%
                        filter(proportion.sd >= threshold)

#total sites after hypo- hyper- filtering
ggplot(proportion_basic_stat, aes(x = proportion.ave)) +
  geom_histogram(binwidth = 0.05, fill = "#009193", color = "black", alpha = 0.7) +
  labs(title = "Average Methylation Proportion across Individuals After Hypo/Hyper/LowV Filter",
       x = "Average Methylation Proportion",
       y = "Frequency") +
  theme_minimal(base_size = 12) +  # Set a base size for text
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
##before sd
ggplot(proportion_basic_stat, aes(x = proportion.sd)) +
  geom_histogram(binwidth = 0.05, fill = "#009193", color = "black", alpha = 0.7) +
  labs(title = "Standard Deviation across Individuals After Hypo/Hyper/LowV Filter",
       x = "Standard Deviation",
       y = "Frequency") +
  theme_minimal(base_size = 12) +  # Set a base size for text
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
```

Added April 2023, I need to remove loci that have one outlier as these warp distribution of methylation for the Matrix QTL analyses (likely other analyses as well, so I thought to remove them and see how the other analyses)

##### Remove Outlier CpG sites

In running some tests for analysis, we realized that there are many sites where one or maybe two individuals have a different methylation value from all other samples. This seems to become correlated to genotyping errors giving many false positives in the analysis. Therefore, I want to make a script that removes these outlier CpG sites (and I want to see the proportion of sites which show this outlier). Therefore, first I want to filter for such outlier CpG sites and then I want to make i) a dataframe which removes these sites. The difference in methylation may be true or due to technical errors,however, without techncial replicates we have no way of knowing.

```r
# Create the initial dataframe
df <- as.data.frame(t(proportion_basic_stat[2:(N+1)]))
colnames(df) <- proportion_basic_stat$site.id

# Function to identify outliers and count them
identify_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 3.5 * IQR
  upper_bound <- Q3 + 3.5 * IQR
  (x < lower_bound) | (x > upper_bound) # Returns logical vector of outliers
}

#Remove these outliers
proportion_basic_stat_rmout <- df %>%
  mutate(across(everything(), ~ifelse(identify_outliers(.x), NA, .x)))

proportion_basic_stat_rmout <- t(proportion_basic_stat_rmout)
colnames(proportion_basic_stat_rmout) <- colnames(proportion_basic_stat[1:86])
#see which were filtered out
# To examine which values were replaced:
outliers_identified <- df %>%
  mutate(across(everything(), identify_outliers))
```

```r
write.table(proportion_basic_stat, file = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(proportion_basic_stat_rmout, file = "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```

Next, we can look at summary statistics using the 3_summary_stat.r script:

```r
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
setwd("/dss/dsshome1/lxc0B/ra52qed/output/1.RRBS/4.0.methylunite/168.0miss/")
#sample names
sample_names <- read.table("/dss/dsshome1/lxc0B/ra52qed/scripts/0.sample.lists/sample.names.168.txt")
#get pop info
info <- sample_names %>% separate(V1, c("pop","year", "subspecies", "id", "bld", "adl", "sex"), "_") 
##Dataset1: No missing data across the entire dataframe
#However, I can still look at average coverage, etc.
load("count.168.0miss.RData") 
#get only total counts
total.counts <- count.168.0miss %>%
 dplyr::select(starts_with("TO_"))
total.counts <- as.data.frame(t(total.counts))
#total.counts <- cbind(sample_names, info$pop, info$subspecies ,total.counts)
#get average coverage per sample
q = c(.25, .5, .75) 
summary.cov.0 <- total.counts %>%
 rowwise %>%
 dplyr::summarize(mean=mean(c_across("V1":"V173933"), na.rm = T),
 min=min(c_across("V1":"V173933"), na.rm = T), 
max=max(c_across("V1":"V173933"), na.rm = T), 
sd=sd(c_across("V1":"V173933"), na.rm = T),
 quant25 = quantile(c_across("V1":"V698368"), na.rm=T, probs = q[1]), 
quant50 = quantile(c_across("V1":"V698368"), na.rm = T, probs = q[2]),
 quant75 = quantile(c_across("V1":"V698368"), na.rm = T, probs = q[3]))
summary.cov.0 <- cbind(sample_names, info$pop, info$subspecies, summary.cov.0) 
save(summary.cov.0, file = "summary.cov.0miss.RData")
```

### 7.1 CpG island dataset

Additionally, I make a dataset of cpg islands. For this I need the positions of the methylation loci, which I grab from the text files as following:

```bash
 awk '{print $1}' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.txt > prop.85.merged.20miss.rmlowv.rmout.positions.tmp

 sed -i 's/_/\t/g' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.positions.tmp
 #deleter site.id header

 awk '{print $1 "\t" $2 "\t" ($2 +1)}' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.positions.tmp > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.positions.bed
```

Then run /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgisland_cpgloci_overlap/1_bedtools.overlap.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=20:00
#SBATCH -J intersect

#run bedtools to find the CpG sites which overlap with repeat regions and CpG islands
#this will later be combined with homer output to look at genomic architecture of CpG sites
bedtools intersect -wa -wb -nonamecheck \
                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.positions.bed \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgisland_cpgloci_overlap/cpgisland_cpgloci_output.85.20miss.bed

#Describing options used in the script
<<COMMENT
  -wa keeps the entries from the A input file
  -wb keeps the entries from the B input file, only for overlapping regions
  -nonamecheck For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. "chr1" vs "chr01".  
  -a input file A (the file used to find overlaps in B)
  -b the file used for finding any overlaps from A
COMMENT
```

which finds overlaps between the identified CpG islands and the called methylation loci. I also want to identify shores (0-2kb from CpG islands) and shelves (2-4kb from CpG islands).

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=20:00
#SBATCH -J intersect

#identify 200bp around the cpgisland...this will define the shore
bedtools slop -b 2000 -i /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged.bed \
                 -g /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/bedtools.cpgi.with.shores.bed

bedtools subtract -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/bedtools.cpgi.with.shores.bed \
                 -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged.bed \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.bed

bedtools intersect -wa -wb -nonamecheck \
                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.positions.bed \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.85.20miss.bed

#identify 2000-4000bp around the cpgisland...this will define the shelf
bedtools slop -b 4000 -i /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged.bed \
                 -g /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/0_rename_chr/chromosomes.genome \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/bedtools.cpgi.with.shelves.bed

bedtools subtract -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/bedtools.cpgi.with.shelves.bed \
                 -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/1_CpGIsland/2_merge/HR_CGI_merged.bed \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.shelves.bed

bedtools subtract -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.shelves.bed \
                 -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.bed \
                 > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shelves.bed

bedtools intersect -wa -wb -nonamecheck \
                  -a /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.positions.bed \
                  -b /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shelves.bed \
                  > /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shelves.85.20miss.bed

#Describing options used in the script
<<COMMENT
slop
    -b the number of basepairs to keep on either side of the genomic regi-
    -i the file to extend 
    -g the genome file
subtract
    -a file to remove regions from
    -b file of regions to remove
intersect
  -wa keeps the entries from the A input file
  -wb keeps the entries from the B input file, only for overlapping regions
  -nonamecheck For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. "chr1" vs "chr01".  
  -a input file A (the file used to find overlaps in B)
  -b the file used for finding any overlaps from A
COMMENT
```

Then run ./0_Pre-processing/4_methyl_calling/1_CpG_islands/0_homer_cpgi_join.Rmd to join the Cpgisland, shores and shelves to the called methylation sites and join this also with homer annotation to have one full document with every. Here, I create one loci based table and one region based table.

````r
---
title: "Defining loci in CpGislands"
author: "Sarah Mueller"
date: "2023-07-20"
output: html_document
---

```{r setup, include=FALSE}
### Libraries
.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))

library(vegan)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(broom)
library(readxl)
library(ggrepel)
library(ggpubr)
library(purrr)
library(wesanderson)
library(readr)
````

## Defining CpG islands from our CpG loci

There are several ways to look at CpG loci and define summary statistics. One option is to look at CpG islands as one unit instead of individual CpG loci. Here, we can look at how intercorrelated the CpG loci are to each other within CpG islands and also what the average methylation across CpG islands is. We can then use the scripts for downstream analysis (DMR detection, PST, PCA,RDA) on the CpG islands and see if this provides us with different information than what we see using only CpG loci. So, the end goal here is to have a dataset that looks similar to prop.168.0miss, but instead of site.id we have a cpgi.id (cpgisland.id).

## Creating CpGi dataset

First, in order to get a file that includes all the called CpG loci and the corresponding CpG islands I used bedtools, see script (/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/1.RRBS/4.4.genomic.regions/0_cpgisland_cpgloci_overlap/1_bedtools.overlap.sh)

##Create CpGi dataset

Now I want to use this file to create a CpGi id so we can later erge all CpG loci that fall in one CpG island.

```{r
#load bedtools output
cpgi.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgisland_cpgloci_overlap/cpgisland_cpgloci_output.85.20miss.bed", header=F)

#create column names
colnames(cpgi.20miss) <- c("chrom.cpgl", "pos.cpgl", "pos1.cpgl", "chrom", "start", "end") 

#merge the chrom, start and end columns into 1 cpgi id column
cpgi.20miss <- cpgi.20miss %>%
                      tidyr::unite(cpgi.id, chrom, start, end, sep = '_') %>%
                      tidyr::unite(site.id, chrom.cpgl, pos.cpgl, sep = '_') 

cpgi.20miss <- cpgi.20miss %>%
                      dplyr::select("site.id", "cpgi.id")
```

Now I want to also merge the shores and shelves data, to have one dataframe with all the different regions included

```{r,
#SHORES
#load bedtools output
cpg.shore.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shores.85.20miss.bed", header=F)

#create column names
colnames(cpg.shore.20miss) <- c("chrom.cpgl", "pos.cpgl", "pos1.cpgl", "chrom", "start", "end") 

#merge the chrom, start and end columns into 1 cpgi id column
cpg.shore.20miss <- cpg.shore.20miss %>%
                      tidyr::unite(cpg.shore.id, chrom, start, end, sep = '_') %>%
                      tidyr::unite(site.id, chrom.cpgl, pos.cpgl, sep = '_') 

cpg.shore.20miss <- cpg.shore.20miss %>%
                      dplyr::select("site.id", "cpg.shore.id")

#SHELVES
#load bedtools output
cpg.shelf.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/2_bedtools/0_cpgi_shores_shelves_overlap/cpgi.shelves.85.20miss.bed", header=F)

#create column names
colnames(cpg.shelf.20miss) <- c("chrom.cpgl", "pos.cpgl", "pos1.cpgl", "chrom", "start", "end") 

#merge the chrom, start and end columns into 1 cpgi id column
cpg.shelf.20miss <- cpg.shelf.20miss %>%
                      tidyr::unite(cpg.shelf.id, chrom, start, end, sep = '_') %>%
                      tidyr::unite(site.id, chrom.cpgl, pos.cpgl, sep = '_') 

cpg.shelf.20miss <- cpg.shelf.20miss %>%
                      dplyr::select("site.id", "cpg.shelf.id")
```

Now I want to merge these cpgi, shores and shelves with the prop.168.20miss dataset and look into how intercorrelated the results are between loci within cpg islands,and how many loci lie outside of cpg islands.

```{r
#load sample names
sample_names <- read.table(file = "/dss/dsshome1/lxc0B/ra52qed/gitlab/swallow.projects/scripts/0_sample.lists/sample.names.85.txt")

#load methylation data
prop.20miss <- read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/0_all_sites/85.20miss/prop.85.merged.20miss.rmlowv.rmout.txt", header=T)

##load homer input 
homer.output.20miss <- read.delim("/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_genome_files/3_CpG_genomic_regions/1_homer/homer.output.85.merged.20miss.rmlowv.rmout.txt")
colnames(homer.output.20miss)[1] <- "PeakID"

#join the 2 dataframes keeping entries from both tables
cpgi.prop.20miss <- full_join(cpgi.20miss, prop.20miss, by = "site.id")
cpgi.prop.20miss <- full_join(cpg.shore.20miss, cpgi.prop.20miss, by = "site.id")
cpgi.prop.20miss <- full_join(cpg.shelf.20miss, cpgi.prop.20miss, by = "site.id")
homer.cpgi.prop.20miss <- left_join(cpgi.prop.20miss, homer.output.20miss, by = c("site.id" = "PeakID"))

homer.cpgi.prop.20miss <- homer.cpgi.prop.20miss[rowSums(is.na(homer.cpgi.prop.20miss[, 5:90])) != ncol(homer.cpgi.prop.20miss[, 5:90]), ]

homer.cpgi.prop.20miss <- homer.cpgi.prop.20miss %>%
  separate(Annotation, into = c("Annotation", "extra"), sep = "\\(", remove = FALSE)

homer.cpgi.prop.20miss <- homer.cpgi.prop.20miss %>%
  mutate(category = case_when(
    !is.na(cpg.shelf.id) ~ "shelf",
    !is.na(cpg.shore.id) ~ "shore",
    !is.na(cpgi.id) ~ "island",
    TRUE ~ "opensea"
  ))

homer.cpgi.prop.20miss <- homer.cpgi.prop.20miss %>%
  mutate(category = paste(category, Annotation, sep = "_")) %>%
  mutate(loci.id = site.id)


homer.cpgi.prop.20miss <- homer.cpgi.prop.20miss %>%
  mutate(region.id = coalesce(cpg.shelf.id, cpg.shore.id, cpgi.id, loci.id))
```

Ok, so now I want to collapse all CpG islands , shores and shelves into one average value per individual per region.

```{r
homer_subset <- homer.cpgi.prop.20miss %>%
        select(c(,91:111))

cpgi_values <- homer.cpgi.prop.20miss %>%
  group_by(cpgi.id) %>%
  dplyr::summarise(across(contains("_ADL_"), mean, na.rm=TRUE))

cpgi_values <- cpgi_values %>%
              left_join(homer_subset, by=c("cpgi.id" = "region.id"))

# Keep only unique entries based on the column 'region.id'
cpgi_values <- cpgi_values %>% rename(region.id = cpgi.id)
cpgi_values <- cpgi_values %>% distinct(region.id, .keep_all = TRUE)


shore_values <- homer.cpgi.prop.20miss %>%
  group_by(cpg.shore.id) %>%
  dplyr::summarise(across(contains("_ADL_"), mean,  na.rm=TRUE))

shore_values <- shore_values %>%
              left_join(homer_subset, by=c("cpg.shore.id" = "region.id"))

# Keep only unique entries based on the column 'region.id'
shore_values <- shore_values %>% rename(region.id = cpg.shore.id)
shore_values <- shore_values %>% distinct(region.id, .keep_all = TRUE)

shelf_values <- homer.cpgi.prop.20miss %>%
  group_by(cpg.shelf.id) %>%
  dplyr::summarise(across(contains("_ADL_"), mean,  na.rm=TRUE))

shelf_values <- shelf_values %>%
              left_join(homer_subset, by=c("cpg.shelf.id" = "region.id"))

# Keep only unique entries based on the column 'region.id'
shelf_values <- shelf_values %>% rename(region.id = cpg.shelf.id)
shelf_values <- shelf_values %>% distinct(region.id, .keep_all = TRUE)

open_sea_values <- homer.cpgi.prop.20miss[,1:111] %>%
              filter(is.na(cpgi.id)) %>%
              filter(is.na(cpg.shore.id)) %>%
              filter(is.na(cpg.shelf.id)) 

# Keep only unique entries based on the column 'region.id'
open_sea_values <- open_sea_values %>% distinct(region.id, .keep_all = TRUE)
open_sea_values <- open_sea_values %>% dplyr::select(-"region.id", -"Chr")
colnames(open_sea_values)[1] <- "region.id"

open_sea_values <- open_sea_values %>%
  dplyr::select(-2:-4)

cat.prop.20miss <- rbind(cpgi_values, shore_values, shelf_values, open_sea_values)

write_delim(cat.prop.20miss,file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/1_CpG_islands/85.20miss/homer.cpgi.prop.85.merged.20miss.rmlowv.rmout.txt", delim="\t", quote_escape = "double")  

homer.cpgi.prop.20miss <- homer.cpgi.prop.20miss[,-c(2:4)]

#There are some places where the shores/shelves overlap because the CpG islands are so close to eachother. Therefore, we want to only keep one shore in this case. Therefore, in the loci_ids we will only keep the assignment to one CpG shore or shelf
homer.cpgi.prop.20miss <- homer.cpgi.prop.20miss %>% distinct(loci.id, .keep_all = TRUE)

# For homer.cpgi.prop.168.20miss
write_delim(homer.cpgi.prop.20miss, file="/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/4_methyl_calling/1_CpG_islands/85.20miss/homer.loci.prop.85.merged.20miss.rmlowv.rmout.txt", delim="\t", quote_escape = "double")  
```

And that concludes pre-prcoessing for the methylation project. See other files for analysis specific methods.