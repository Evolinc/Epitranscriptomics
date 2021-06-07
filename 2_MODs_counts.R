#!/usr/bin/Rscript
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#load basic packages
library(getopt)
library(tidyr)
library(reshape2)
library(dplyr)


## -----------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'ratio', 'r', 1, "character",
    'output1', 'o1', 1, "character",
    'output2', 'o2', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)


##help info---------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$ratio) || is.null(args$output1) || is.null(args$output2)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}



## -----------------------------------------------------------------------------------------------------------------------------------
#Load files
Read_Meta <- read.csv(args$ratio, header = T)

## ---------------------------------------calculate mods numbers (all and sep) under same gene-----------------------------------------

#calculate respective mods numbers under same gene
MODs_count_sep <- Read_Meta[,c(4,6,7)]%>%
  group_by(Gene.ID,MODs_type,Genomic.elemnet) %>%
  dplyr::mutate(Counts =  n())

#Re-structure mods
MODs_count_sep_dcast <- dcast(MODs_count_sep, Genomic.elemnet + Gene.ID ~ MODs_type)


#calculate all mods numbers under same gene
MODs_count_all <- Read_Meta[,c(4,6,7)]%>%
  group_by(Gene.ID,Genomic.elemnet) %>%
  dplyr::mutate(Counts =  n())

MODs_count_all <- unique(MODs_count_all[,c(1,2,4)])

## -----------------------------------------------------------------------------------------------------------------------------------
#save results
write.csv(MODs_count_sep_dcast, args$output1, row.names = F) 
write.csv(MODs_count_all, args$output2, row.names = F) 



