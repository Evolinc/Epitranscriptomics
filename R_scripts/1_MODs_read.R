#!/usr/bin/Rscript

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(getopt)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'reads', 'r', 1, "character",
    'mods', 'm', 1, "character",
    'output_table', 't', 1, "character",
    'output_meta', 'o', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)


##help info-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$reads) || is.null(args$mods) || is.null(args$output_table) || is.null(args$output_meta)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load mapping reads file
MODs_mapping <- read.table(args$reads, header = T)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load mods annotation file
MODs_anno <- read.delim(args$mods, "\t",header = FALSE)

#Reformat the Mods annotaion
MODs_anno_clean <- MODs_anno[,c(2,3,5,8,10,11,24)]
names(MODs_anno_clean)[1:7] <- c("Chr", "bp", "Strand", "Region", "Distance_TSS", "Gene.ID","MODs_type")
MODs_anno_clean[,c(2,5)] <- sapply(MODs_anno_clean[,c(2,5)],as.integer)
MODs_anno_clean <- MODs_anno_clean %>% separate(Region, c("Genomic.elemnet"), sep =" ")
MODs_anno_clean <- MODs_anno_clean %>% separate(Gene.ID, "Gene.ID", sep =":")

## ----------------------------------------------------------merge ratio and loci information to certain gene annotations-----------------------------------------------------------------
#Join dataset (Annotation and SNP table)
##non-genic part will be removed
Read_meta <- subset(full_join(MODs_anno_clean ,MODs_mapping ,by = "bp")[,c(1,2,3,4,5,6,7,10,11,12,13,14,15,16,17,18,19,20)], Gene.ID != "NA")

#Calculate ratio for each site
Read_meta$ratio <- Read_meta$nonref/(Read_meta$nonref + Read_meta$ref)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(Read_meta, args$output_table, row.names=FALSE) 



## -------------------------------------------------------------Calculate modififed ratio for each gene----------------------------------------------------------------------------------
Read_meta_gene_ratio <- Read_meta[,c(4,6,7,13,14)]%>%
group_by(Gene.ID,MODs_type,Genomic.elemnet) %>%
dplyr::mutate(nonref_sum =  sum(nonref))
    
Read_meta_gene_ratio <- Read_meta_gene_ratio%>%
group_by(Gene.ID,MODs_type,Genomic.elemnet) %>%
dplyr::mutate(ref_sum =  sum(ref))

#calculate the ratio of reads under certain genes
Read_meta_gene_ratio$Gene_ratio <- Read_meta_gene_ratio$nonref_sum/(Read_meta_gene_ratio$nonref_sum + Read_meta_gene_ratio$ref_sum)


#Re-structure mods using Dcast function
Read_meta_gene_ratio_dcast <- dcast(unique(Read_meta_gene_ratio[,c(1,2,3,8)]), Genomic.elemnet+Gene.ID~MODs_type)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(Read_meta_gene_ratio_dcast, args$output_meta, row.names=FALSE) 



