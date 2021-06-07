#!/usr/bin/Rscript

## -----------------------------------------------------------------------------------------------------------------------------------
library(plyr)
library(dplyr)
library(getopt)

## -----------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'mods_count', 'mc', 1, "character",
    'mods_ratio', 'mr', 1, "character",
    'synteny_table', 't', 1, "character",
    'Synteny_count', 'c', 1, "character",
    'Synteny_ratio', 'r', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)


##help info---------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$mods_count) || is.null(args$mods_ratio) || is.null(args$synteny_table) || is.null(args$Synteny_count) || is.null(args$Synteny_ratio)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

## -----------------------------------------------------------------------------------------------------------------------------------

Plot_path <- '/home/liangyu/1_Project/3_RNA-mods/3_Plot/3_synteny'

## -----------------------------------------------------------------------------------------------------------------------------------
#Load modification summary sheet
MODs_count <- read.csv(args$mods_count, header = T)
MODs_count[is.na(MODs_count)] <- 0
#merge three genomic elements together
MODs_count <- ddply(MODs_count, "Gene.ID", numcolwise(sum))


MODs_ratio <- read.csv(args$mods_ratio, header = T)
MODs_ratio[is.na(MODs_ratio)] <- 0
MODs_ratio <- ddply(MODs_ratio, "Gene.ID", numcolwise(sum))


#Load syntenic genes list
synteny.df <- read.csv(args$synteny_table, header = F)
Synteny_1 <- data.frame("Gene.ID" = synteny.df[,1])
Synteny_2 <- data.frame("Gene.ID" = synteny.df[,2])


## -----------------------------Join two sets of the synteny table to modification count sheet----------------------------------------
Synteny_count_1 <- left_join(Synteny_1, MODs_count, by = "Gene.ID")
Synteny_count_2 <- left_join(Synteny_2, MODs_count, by = "Gene.ID")

Synteny_count <- cbind(Synteny_count_1,Synteny_count_2)


## -----------------------------Join two sets of the synteny table to modification ratio sheet----------------------------------------
Synteny_ratio_1 <- left_join(Synteny_1, MODs_ratio, by = "Gene.ID")
Synteny_ratio_2 <- left_join(Synteny_2, MODs_ratio, by = "Gene.ID")

Synteny_ratio <- cbind(Synteny_ratio_1,Synteny_ratio_2)

## -----------------------------------------------------------------------------------------------------------------------------------
write.csv(Synteny_count, args$Synteny_count, row.names = F) 
write.csv(Synteny_ratio, args$Synteny_ratio, row.names = F) 
