## HOMER Annotation

I wanted to identify in what functional categories the CpG sites I identified came from.

I decided to use [HOMER]([Homer Software and Data Download](http://homer.ucsd.edu/homer/motif/)) to do this

First, homer needs an input file for a custom .fasta file (one that is not in the reference bank). Additionally, I need a bed file that has all of my CpG sites...it looks like this:

chr start end id strand
chr01 8675 8675 chr01_8675 1
chr01 10873 10873 chr01_10873 1
chr01 10884 10884 chr01_10884 1
chr01 10891 10891 chr01_10891 1
chr01 10977 10977 chr01_10977 1

and is created using 0_homer_input.sh

I then run homer using 1_homer.sh. I define my reference genome and annotate peaks (in my case CpG sites).

Visualization is done using the R scripts within this, but these are not polished scripts, rather just some unfinished plots, etc.