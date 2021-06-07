#!/usr/bin/Rscript

## -----------------------------------------------------------------------------------------------------------------------------------
library(getopt)
library(biomaRt)
library(org.At.tair.db)
library(clusterProfiler)
library(dplyr)

##Plot packages
library(DOSE)
library(ggplot2)
library(ggnewscale)


## -----------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'gene_list', 'g', 1, "character",
    'GO_table', 'gt', 1, "character",
    'KEGG_table', 'kt', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)


##help info---------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$gene_list) || is.null(args$GO_table) || is.null(args$KEGG_table)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

## -----------------------------------------------------------------------------------------------------------------------------------
Plot_path <- '/home/liangyu/1_Project/3_RNA-mods/3_Plot/4_GO'


## ---------------------------------------Load Genome annotations--------------------------------------------------------
ensembl <- useEnsembl(biomart = "genes")
ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")
ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart", 
                                         dataset = "athaliana_eg_gene")

## -----------------------------------------------------------------------------------------------------------------------------------
#Load the annotation of MODs associated with gene IDs
MODs_list <- read.csv(args$gene_list, header = T)


## ----------------------------------------Check types of modifications within MODs_list-----------------------------------------------
temp <- unique(MODs_list$MODs_type)

#Build dataframe for all MODs appeared in list
MODs.df <- data.frame(name = paste0(temp,"_GO"),
            Type = temp)

## -------------------------------------------GO Enrichment for listed modified genes----------------------------------------------------


for (i in 1:nrow(MODs.df)){
   df <- MODs_list[MODs_list$MODs == MODs.df[i,2],]
  
   #Transform the Gene ID to entrezgene_id
   df_ID = paste0(MODs.df[i,2],"_ID")
   assign(df_ID, getBM(attributes=c("entrezgene_id"),filters = "ensembl_gene_id",values = df[,1], mart = ensembl_arabidopsis))
   
   ##Read (get) the df_ID dataframe from global Env (IMPORTANT)
   df_ID <- get(paste0(MODs.df[i,2],"_ID"))
   #Enrichement for biological process
   df_BP = paste0(MODs.df[i,2], "_BP")
   assign(df_BP, enrichGO(gene = df_ID[,1], OrgDb = org.At.tair.db, keyType = "ENTREZID",
                       ont = "BP",qvalueCutoff = 0.2, pvalueCutoff = 0.05))
   
   #Enrichement for Cellular comonent
   df_CC = paste0(MODs.df[i,2], "_CC")
   assign(df_CC, enrichGO(gene = df_ID[,1], OrgDb = org.At.tair.db, keyType = "ENTREZID",
                       ont = "CC",qvalueCutoff = 0.2, pvalueCutoff = 0.05))
   
   #Enrichement for Molecular Function
   df_MF = paste0(MODs.df[i,2], "_MF")
   assign(df_MF, enrichGO(gene = df_ID[,1],OrgDb = org.At.tair.db, keyType = "ENTREZID",
                       ont = "MF",qvalueCutoff = 0.2, pvalueCutoff = 0.05))

   #Assign variable for looping
   assign(MODs.df[i,1],df)
}


## ----------------------------------------Export the summary associate genes from Molecullar Function-------------------------------------------

#Create list of MF ontology list
rm(df_MF)
MF_list <- ls(pattern = "_MF")

#Transform each categories into dataframe
for (i in 1:length(MF_list)){
   if(class(get(MF_list[i])) == "enrichResult"){
      temp.df.summary <- as.data.frame(setReadable(get(MF_list[i]), 'org.At.tair.db', 'ENTREZID'))
      if (nrow(temp.df.summary)!=0){
         temp.df.summary$MODs_type <- as.character(paste0(gsub('.{3}$', '', MF_list[i]))) 
         assign(paste0(MF_list[i],"_summary"),temp.df.summary)
      }
   }
}
 
rm(temp.df.summary)
#Combine all GO_MF 
GO_MF <- bind_rows(lapply(ls(pattern = "_MF_summary"), get))
GO_MF$GO_type <- "Molecular Function"


## ----------------------------------------Export the summary associate genes from Biological Process-------------------------------------------

#Create list of BP ontology list
rm(df_BP)
BP_list <- ls(pattern = "_BP")

#Transform each categories into dataframe
for (i in 1:length(BP_list)){
   if(class(get(BP_list[i])) == "enrichResult"){
      temp.df.summary <- as.data.frame(setReadable(get(BP_list[i]), 'org.At.tair.db', 'ENTREZID'))
      if (nrow(temp.df.summary)!=0){
         temp.df.summary$MODs_type <- as.character(paste0(gsub('.{3}$', '', BP_list[i]))) 
         assign(paste0(BP_list[i],"_summary"),temp.df.summary)
      }
   }
}

rm(temp.df.summary)
#Combine all GO_MF 
GO_BP <- bind_rows(lapply(ls(pattern = "_BP_summary"), get))
GO_BP$GO_type <- "Biological Process"



## ----------------------------------------Export the summary associate genes from Cellular Component-------------------------------------------

#Create list of CC ontology list
rm(df_CC)
CC_list <- ls(pattern = "_CC")

#Transform each categories into dataframe
for (i in 1:length(CC_list)){
  if(class(get(CC_list[i])) == "enrichResult"){
    temp.df.summary <- as.data.frame(setReadable(get(CC_list[i]), 'org.At.tair.db', 'ENTREZID'))
    if (nrow(temp.df.summary)!=0){
      temp.df.summary$MODs_type <- as.character(paste0(gsub('.{3}$', '', CC_list[i]))) 
      assign(paste0(CC_list[i],"_summary"),temp.df.summary)
    }
  }
}

rm(temp.df.summary)
#Combine all GO_MF 
GO_CC <- bind_rows(lapply(ls(pattern = "_CC_summary"), get))
GO_CC$GO_type <- "Cellular Component"


## ----------------------------------------------Combine enriched GOs for each type of modifications--------------------------------------------------

GO_combine <- bind_rows(GO_CC,GO_BP,GO_MF)
write.csv <- (GO_combine,args$GO_table, row.names=FLASE)



## -------------------------------------------KEGG Enrichment for listed modified genes----------------------------------------------------

KEGG.df <- data.frame(name = paste0(temp,"_gene"),
            Type = paste0(temp,"_KEGG"))

#Perform KEGG analylysis for each gene class
for (x in 1:nrow(KEGG.df)){
    #Load Gene list for each type of the mods
    kk <- enrichKEGG(gene = get(KEGG.df[x,1])[,1],
                 organism = 'ath', 
                 pvalueCutoff = 0.05)
    
    Kk_summary <- as.data.frame(kk)
    assign(KEGG.df[x,2],kk)
    assign(paste0(KEGG.df[x,2], "_summary"),Kk_summary)
}


#list modifications along with summary output
rm(Kk_summary)
sum_list <- ls(pattern = "summary")


#Add additional column to specify modification types
for (i in 1:length(sum_list)){
   temp.df <- get(sum_list[i])
   if (nrow(temp.df)!=0) {
   temp.df$MODs_type<- as.character(paste0(gsub('.{18}$', '', sum_list[i]))) 
   assign(sum_list[i],temp.df)
   }
}

## ----------------------------------------------Combine enriched KEGGs for each type of modifications--------------------------------------------------

#Build a list of dataframe contains enrichment
merge_list <- lapply(ls(pattern = "_summary"), get)
KEGG_combine <- bind_rows(merge_list)

write.csv <- (KEGG_combine,args$KEGG_table, row.names=FLASE)


