#!/bin/bash

# Script to process RNA modificaiton data for single sample
# Multiple datasets will be generated under the same output directory

usage() {
      echo ""
      echo "Usage : sh $0 -b MODs_bed_file -r reads_depth_file -s synteny -g gene_gff -G ref_genome -e Repeat_gff -m motif -o output_dir -h"
      echo ""

cat <<'EOF'

  -b </path/to/MODs bed output file>

  -r </path/to/MODs reads depth file>

  -s </path/to/ gene synteny file>

  -g </path/to/ gene gff file>

  -G </path/to/ reference genome fasta file>

  -e </path/to/Repeat element gff file>

  -m <maximum numbers of motifs for MEME enrichment>

  -o </path/to/output>

  -h Show this usage information

EOF
    exit 0
}

while getopts ":b:r:s:g:G:e:m:o:h:" opt; do
  case $opt in
    b)
     MODs_bed_file=$OPTARG
      ;;
    r)
     reads_depth_file=$OPTARG
      ;;
    s)
     synteny=$OPTARG
      ;;
    g)
     gene_gff=$OPTARG 
      ;;  
    G)
     ref_genome=$OPTARG 
      ;;   
    e)
     Repeat_gff=$OPTARG 
      ;;
    m)
     motif=$OPTARG
      ;;
    o)
     output_dir=$OPTARG
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
echo "******************STEP 1 Pre-processing modifications IDs derived from HAMR......******************"

#if [ ! -e "$ref_genome.dict" ]; then
  #java -jar /home/liangyu/bin/picard.jar CreateSequenceDictionary -R $ref_genome -O ${output_dir}/$ref_genome.dict
#else

bed_file="$(echo ${MODs_bed_file} |sed -e 's@.*/@@')" ;
cp $MODs_bed_file $bed_file 
echo "The bed file been processed is $bed_file ......"

#Split two strians of mods for annotation
grep "+" $bed_file > ${output_dir}/$bed_file.plus.txt
grep "-" $bed_file > ${output_dir}/$bed_file.minus.txt

#run hommer annotate peak to annotate coordinates of mods
perl annotatePeaks.pl ${output_dir}/$bed_file.plus.txt $ref_genome -gff $gene_gff -strand + -annStats ${output_dir}/$bed_file.plus_stats > ${output_dir}/$bed_file.plus_anno.txt
perl annotatePeaks.pl ${output_dir}/$bed_file.minus.txt $ref_genome -gff $gene_gff -strand - -annStats ${output_dir}/$bed_file.minus_stats > ${output_dir}/$bed_file.minus_anno.txt

#subtract the SRR ID
SRR=$(echo ${bed_file} | cut -d '.' -f 1)

#merge the plus and minus strand sequences
cat  ${output_dir}/$bed_file.plus_anno.txt ${output_dir}/$bed_file.minus_anno.txt | sort -k1,2 -k2,3 | grep "Chr" > ${output_dir}/$SRR.merge_anno.txt

#Reformat the BED file
cat  $bed_file | sort -k1,1 -k2,2 | grep "Chr" > ${output_dir}/$SRR.bed

#merge (paste) the BED file and annotation file
paste ${output_dir}/$SRR.merge_anno.txt ${output_dir}/$SRR.bed > ${output_dir}/$SRR.anno

#Reformat the final annotation clean file
sed -i 's/ID=//g' ${output_dir}/$SRR.anno

#clean the final annotation file
grep -v 'cmd=annotatePeaks.pl' ${output_dir}/$SRR.anno | cut -f2,3,5,10,11,23,24 > ${output_dir}/$SRR.anno.temp.txt
sed -i 's/.exon/:exon/g' ${output_dir}/$SRR.anno.temp.txt

grep -v 'cmd=annotatePeaks.pl' ${output_dir}/$SRR.anno | cut -f8 | cut -d ' ' -f1 > ${output_dir}/$SRR.anno.element.txt
paste ${output_dir}/$SRR.anno.temp.txt ${output_dir}/$SRR.anno.element.txt > ${output_dir}/$SRR.anno.clean

#Rename modification type to prevent parsing error 
sed -i 's/|/./g' ${output_dir}/$SRR.anno.clean
sed -i 's/|/./g' $reads_depth_file

rm ${output_dir}/*.txt
rm ${output_dir}/*.bed

echo "******************STEP 1 Annotations of modifications completed******************"
ELAPSED_TIME_1=$(($SECONDS - $START_TIME_1))
echo "Elapsed time for STEP 1 is" $ELAPSED_TIME_1 "seconds"



START_TIME_2=$SECONDS
echo "******************STEP 2 Merging mapped reads and ratio of modified reads ......******************"

Rscript 1_MODs_read.R --reads $reads_depth_file --mods ${output_dir}/$SRR.anno.clean --output_table ${output_dir}/$SRR.reads.ratio.csv --output_meta ${output_dir}/$SRR.reads.stats.csv

ELAPSED_TIME_2=$(($SECONDS - $START_TIME_2))
echo "Elapsed time for STEP 2 is" $ELAPSED_TIME_2 "seconds" 

START_TIME_3=$SECONDS
echo "******************STEP 3 Generating modification counts per gene under each modification type ......******************"

Rscript 2_MODs_counts.R --ratio ${output_dir}/$SRR.reads.ratio.csv --output1 ${output_dir}/$SRR.counts.sep.csv --output2 ${output_dir}/$SRR.counts.all.csv

ELAPSED_TIME_3=$(($SECONDS - $START_TIME_3))
echo "Elapsed time for STEP 3 is" $ELAPSED_TIME_3 "seconds"

START_TIME_4=$SECONDS
echo "******************STEP 4 Performing profiling of mods for syntenic genes ......******************"

Rscript 3_MODs_synteny.R --mods_count $SRR.counts.sep.csv --synteny_table $synteny --synteny_count $SRR.syntney.counts.csv

ELAPSED_TIME_4=$(($SECONDS - $START_TIME_4))
echo "Elapsed time for STEP 4 is" $ELAPSED_TIME_4 "seconds"

START_TIME_5=$SECONDS
echo "******************STEP 5 Screening ROIs for MEME enrichment (50% highly-modified genes) ......******************"

Rscript 4_MODs_MEME.R --input ${output_dir}/$SRR.reads.ratio.csv --output ${output_dir}/$SRR.MEME.input.csv

ELAPSED_TIME_5=$(($SECONDS - $START_TIME_5))
echo "Elapsed time for step 5 is" $ELAPSED_TIME_5 "seconds"


