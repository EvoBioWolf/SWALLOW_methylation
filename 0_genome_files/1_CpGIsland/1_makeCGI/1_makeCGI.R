getwd()
setwd("/dss/dsshome1/lxc0B/ra52qed/scripts/0.genome.files/1.makeCGI")

#install.packages("makeCGI_1.3.4.tar.gz", repos = NULL, type = "source")

library(makeCGI)

#Set up default parameters:
.CGIoptions=CGIoptions()
#If the inputs are fa files, users need to create a directory called "rawdata" under current working directory and save all fa files there, then speicify the input type and species name. For example,
.CGIoptions$rawdat.type="txt"
.CGIoptions$species="Hirundo_rustica"

 makeCGI(.CGIoptions)
                     
