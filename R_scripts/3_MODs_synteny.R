#!/usr/bin/Rscript

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(getopt)
library(reshape2)
library(dplyr)
library(tidyr)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


command = matrix(c(
    'mods_count', 'm', 1, "character",
    'synteny_table', 't', 1, "character",
    'synteny_count', 'c', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)


##help info------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if (!is.null(args$help) || is.null(args$mods_count) || is.null(args$synteny_table) || is.null(args$synteny_count)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

## ------------------------------------------------------------Load modification meta_table------------------------------------------------------------------------------------

#Load modification meta table
MODs_meta_sep <- read.csv(args$mods_count, header = T)

#combine sub -modification type from different genomic elements but under the same gene
MODs_meta_sep <- MODs_meta_sep %>% group_by(Gene.ID)%>% dplyr::mutate(sum(D))%>% dplyr::mutate(sum(Y))%>%
  dplyr::mutate(sum(i6A.t6A))%>% dplyr::mutate(sum(m1A.m1I.ms2i6A))%>% dplyr::mutate(sum(m1G))%>%
  dplyr::mutate(sum(m2G.m22G))%>% dplyr::mutate(sum(m3C))

MODs_meta_sep <- unique(data.frame(MODs_meta_sep[,2], MODs_meta_sep[,10:16]))
MODs_meta_sep$Total <- rowSums(MODs_meta_sep[,2:8])

## -------------------------------------------------------------Load syntenic gene list generated by CoGe------------------------------------------------------------------------
#Load syntenic genes list
synteny <- read.csv(args$synteny_table, header = F, na.strings=c("","NA"))

## -------------------------------------------------------- Combine and parse the dataframe ----------------------------------------------------------------------------------------


for (i in 1:ncol(synteny)){
  synteny.df <- left_join(data.frame(Gene.ID = synteny[,i]),MODs_meta_sep, by = "Gene.ID")
  assign(paste0("sub_synteny_",i),synteny.df)
}

#merge all sub_synteny file derived from each column 
merge_list <- lapply(ls(pattern = "sub_synteny_"), get)
synteny_combine <- bind_cols(merge_list)

#Replace all na counts as zero
synteny_combine[is.na(synteny_combine)] = 0

#Remove rows without any modifications 

#calculate the maximum colnumbers need to be sum up
select_col <- seq(9,ncol(synteny_combine),9)

synteny_combine$Total <- rowSums(synteny_combine[,c(select_col)])
synteny_clean <- synteny_combine[synteny_combine$Total != 0,]


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

write.csv(synteny_clean,args$synteny_count, row.names = F)


