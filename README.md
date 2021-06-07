# Epitranscriptomics
The High-throughput analysis pipeline of RNA modifications

## Introduction
Post-transcriptional mRNA modifications play substantial roles of regulating biological processes in plants. The workflow here is to effectively classify, characterize, and compare a variety of RNA modifications identified from HAMR workflow derived from RNA-sequencing data. 

## General analysis of each dataset
### Rationale
Each raw dataset will be processed under the same parameter to characterize distribution pattern and other genomic features of different types of modifications.
Five general topics will be performed for each dataset as follow:

1. Genomic annotation of modifications based on gene annotation from respective species genomes.
2. Calculations of numbers of modifications and portions of modified reads over modified reads from both gene and single locus level
3. Comparisons of modification numbers from syntenic genes
4. Identification of enriched motif(s) for each type of modification
5. Comparison of gene density, modification numbers, and adjacent transposable elements frequency from the same dataset

### input data
1. Reference genome sequences (FASTA)
2. Reference genome gene annotation (gff/gff3/gtf)
3. Known modifications position identified by HAMR (BED)
4. Known modifications mapping reads counts and modified reads counts identified by HAMR (txt)
5. Syntenic gene list
6. Transposable elements annotation (gff)

## Integrated analysis for certain experimental setup
### Rationale
General analysis for each dataset will be integrated to address certain biological questions of RNA modifications
The capability of pipelines can be summarized as:

1. Compare differentially modified genes (DMGs) from multiple experiments
2. Generate enriched gene ontology and pathways for certain gene lists
3. Comparisons of modifications of syntenic genes across species


### input data
1. List of modified genes
2. Statistics of modification derived from general analysis
3. Comparative genomics results (gene synteny) from inter-species comparisons


A schematic workflow of integrated analyis was shown:
![MODs_pipeline](https://user-images.githubusercontent.com/69836931/121089630-72904e80-c7b5-11eb-9b87-c3a36dc83bac.png)

