#!/usr/bin/Rscript

## -----------------------------------------------------------------------------------------------------------------------------------
library(plyr)
library(dplyr)
library(getopt)

## -----------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'input', 'i', 1, "character",
    'output', 'o', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)


##help info---------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$input) || is.null(args$output)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

## -----------------------------------------------------------------------------------------------------------------------------------
Read_meta <- read.csv(args$input, header = T )


## -----------------------------------------------------------------------------------------------------------------------------------
#calculate repeat counts for mods of each gene
Read_meta <- Read_meta %>% group_by(MODs_type,Gene.ID) %>%  dplyr::mutate(Counts =n())

#Subtract gene-level mods summary
MODs_gene <- unique(Read_meta[,c(6,7,20)])

#calculate genes been modified each categories
MODs_df <- MODs_gene %>% group_by(MODs_type) %>%  dplyr::mutate(Genes.modified = n())
MODs_df <- unique(MODs_df[c(2,4)])

#calculate the quantile of gene-level mods for each modification type
quantile <- data.frame(do.call("rbind",tapply(MODs_gene$Counts, MODs_gene$MODs_type,quantile)))
quantile <- tibble::rownames_to_column(quantile, "MODs_type")

#Joint genes been modified each categories and numbers of mods each gene
joint_df <- full_join(quantile,MODs_df,"MODs_type")

## -----------------------------------------------------------------------------------------------------------------------------------

#Take the modifications with mods numbers greater than median of counts (can be 75% percentile)
for (i in 1:nrow(joint_df)){
  df <- MODs_gene[MODs_gene$MODs_type == joint_df[i,1] & MODs_gene$Counts > joint_df[i,4],]
  assign(paste0(joint_df[i,1],"_enrich"),df)
}

#subtract the dataframe
Read_meta_join <- Read_meta[,c(6,1,2,7)]

#define the new dataframe for searching selected modifications
search_df <- data.frame("MODs_type" =  paste0(joint_df[,1]),"Enrich_type" = paste0(joint_df[,1],"_enrich"),"Enrich_corrdinate" = paste0(joint_df[,1],"_coordinate"))

for (n in 1:nrow(search_df)){
  coordinate.df <- left_join(get(search_df[n,2])[,1],Read_meta_join[Read_meta_join$MODs_type == search_df[n,1],],"Gene.ID")
  assign(search_df[n,3],coordinate.df)
}

#Combine all coordinate dataframe
MODs_join <- bind_rows(D_coordinate,Y_coordinate,m1A.m1I.ms2i6A_coordinate,
                       m3C_coordinate,m2G.m22G_coordinate,m1G_coordinate,i6A.t6A_coordinate)

MODs_join <- MODs_join %>% group_by(Gene.ID,MODs_type) %>%  dplyr::mutate(Count = n())
MODs_join <- MODs_join %>% group_by(Gene.ID,MODs_type) %>%  dplyr::mutate(min = min(bp))
MODs_join <- MODs_join %>% group_by(Gene.ID,MODs_type) %>%  dplyr::mutate(max = max(bp))
MODs_join <- MODs_join %>% group_by(Gene.ID,MODs_type) %>%  dplyr::mutate(Range = max(bp)-min(bp))

MODs_join_stats <- unique(MODs_join[c(1,2,4,5,6,7,8)])
MODs_join_stats$bp_per_MODs <- MODs_join_stats$Range/MODs_join_stats$Count


#Took the genes associated with a lower bp_per_mods ratio for enrichement
MODs_join_stats_cutoff <- MODs_join_stats[MODs_join_stats$bp_per_MODs < 50,]
#Export files for MEME enrichment
write.table(MODs_join_stats_cutoff, args$output, row.names=FALSE,col.names=FALSE)
