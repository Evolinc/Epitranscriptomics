#!/bin/bash

# Script to generate enrichment of gene ontology (GO) and KEGG pathways 
# A list of genes for enrichment and speceis information are required 

usage() {
      echo ""
      echo "Usage : sh $0 -l gene_list -g GO_table -K KEGG_table"
      echo ""

cat <<'EOF'

  -l </path/to/gene list file>

  -g </path/to/ saving final GO enrichment output>

  -k </path/to/ saving final KEGG enrichment output>

  -h Show this usage information

EOF
    exit 0
}

while getopts ":b:r:s:g:G:e:m:o:h:" opt; do
  case $opt in
    l)
     gene_list=$OPTARG
      ;;
    g)
     GO_output=$OPTARG
      ;;
    k)
     KEGG_output=$OPTARG
      ;;
    h)
     usage
     exit 1
      ;;      
    \?)
     echo "Invalid option: -$OPTARG" >&2
     exit 1
      ;;
    :)
     echo "Option -$OPTARG requires an argument." >&2
     exit 1
      ;;
  esac
done

echo "BEGIN!"
echo `date`

START_TIME=$SECONDS

START_TIME_1=$SECONDS
echo "******************Performing GO and KEGG enrichment......******************"

Rscript $pipe_dir/4_MODs_GO.R --gene_list $gene_list --GO_table $GO_output --KEGG_table $KEGG_output

echo "**************************Enrichemnt compelted*****************************"
ELAPSED_TIME_1=$(($SECONDS - $START_TIME_1))
echo "Elapsed time for GO and KEGG enrichment is" $ELAPSED_TIME_1 "seconds"