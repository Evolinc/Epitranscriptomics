#!/bin/bash

# Script to process RNA modificaiton data for single sample
# Multiple datasets will be generated under the same output directory

usage() {
      echo ""
      echo "Usage : sh $0 -b modsbed -r reads_depth -s synteny -g gene_gff -G genome -e repeat -m motif -o outputdir"
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
     outputdir=$OPTARG
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

java -jar picard.jar CreateSequenceDictionary -R $ref_genome -O $ref_genome.dict

#Split two strians of mods for annotation
grep "+" $MODs_bed_file > $MODs_bed_file.plus
grep "-" $MODs_bed_file > $MODs_bed_file.minus

#run hommer annotate peak to annotate coordinates of mods
perl ./annotatePeaks.pl $MODs_bed_file.plus $ref_genome -gff $gene_gff -strand + -annStats $MODs_bed_file.plus_stats > $MODs_bed_file.plus_anno
perl ./annotatePeaks.pl $MODs_bed_file.minus $ref_genome -gff $gene_gff -strand - -annStats $MODs_bed_file.minus_stats > $MODs_bed_file.minus_anno

#merge annotation file and BED file
cat  $MODs_bed_file.????_anno | sort -k1,2 -k2,3 | grep "Chr" > $MODs_bed_file.merge_anno
cat  ${MODdir}/$sample.* | sort -k1,1 -k2,2 | grep "Chr" > ${outdir}/$sample.merge.bed

paste ${outdir}/$sample.merge.anno.txt ${outdir}/$sample.merge.bed > $outdir/${sample}.anno.clean.txt
sed -i 's/(ID=//g' | sed -i 's/ID=//g' $outdir/${sample}.anno.clean.txt
grep -v 'cmd=annotatePeaks.pl' $outdir/${sample}.anno.clean.txt | cut -f1,2,3,4,5,6,10,11,23,24 > $MODs_bed_file.clean

#Rename modification type to prevent parsing error 
sed -i 's/|/./g' $MODs_bed_file.clean
sed -i 's/|/./g' $reads_depth_file

echo "******************STEP 1 Annotations of modifications completed******************"
ELAPSED_TIME_1=$(($SECONDS - $START_TIME_1))
echo "Elapsed time for STEP 1 is" $ELAPSED_TIME_1 "seconds"

START_TIME_2=$SECONDS
echo "******************STEP 2 Merging mapped reads and ratio of modified reads ......******************"

Rscript ./1_MODs_read.R --reads $reads_depth_file --mods $MODs_bed_file.clean --output_table ${outputdir}/$reads_depth_file.reads.ratio.csv --output_meta ${outputdir}/$reads_depth_file.reads.stats.csv

ELAPSED_TIME_2=$(($SECONDS - $START_TIME_2))
echo "Elapsed time for STEP 2 is" $ELAPSED_TIME_2 "seconds" 

START_TIME_3=$SECONDS
echo "******************STEP 3 Generating modification counts per gene under each modification type ......******************"

Rscript ./2_MODs_counts.R --ratio ${outputdir}/$reads_depth_file.reads.ratio.csv --output1 ${outputdir}/$MODs_bed_file.counts.sep.csv --output2 ${outputdir}/$MODs_bed_file.counts.all.csv

ELAPSED_TIME_3=$(($SECONDS - $START_TIME_3))
echo "Elapsed time for STEP 3 is" $ELAPSED_TIME_3 "seconds"
	

START_TIME_4=$SECONDS
echo "******************STEP 4 Performing profiling of mods for syntenic genes ......******************"

Rscript ./3_MODs_synteny.R --mods_count ${outputdir}/$MODs_bed_file.counts.sep.csv --mods_ratio ${outputdir}/$reads_depth_file.reads.stats.csv --synteny_table $synteny --Synteny_count ${outputdir}/$reads_depth_file.syntney.counts.csv --Synteny_ratio ${outputdir}/$reads_depth_file.syntney.ratio.csv

ELAPSED_TIME_4=$(($SECONDS - $START_TIME_4))
echo "Elapsed time for STEP 4 is" $ELAPSED_TIME_4 "seconds"

#STEP4.Screen ROIs for MEME enrichment (50% highly-modified genes)

START_TIME_5=$SECONDS
echo "******************STEP 5 Screening ROIs for MEME enrichment (50% highly-modified genes) ......******************"

Rscript ./4_MODs_MEME.R --input ${outputdir}/$reads_depth_file.reads.ratio.csv --output ${outputdir}/$reads_depth_file.MEME.input.txt

ELAPSED_TIME_5=$(($SECONDS - $START_TIME_5))
echo "Elapsed time for step 5 is" $ELAPSED_TIME_5 "seconds"


START_TIME_6=$SECONDS
echo "******************STEP 6 Performing motif enrichment by MEME ......******************"

cat ${outputdir}/$reads_depth_file.MEME.input.txt | cut -d ' ' -f 1,2,3,4,5,6 | sed 's/"//g' | sed 's/Chr//g' >  ${outputdir}/$reads_depth_file.MEME.input.txt

for n in i6A.t6A m3C D m1G m1A.m1I.ms2i6A Y m2G.m22G;
do
  cd $outputdir
  mkdir $n
  cd $n
  
  #divide annotations into different classes
  grep $n ${outputdir}/$reads_depth_file.MEME.input.txt > $outputdir/${n}/$reads_depth_file.MEME._$n.txt

  #Generate corrdinate list for samotools
  awk '{ print $2 ":" $5 "-" $6, $1 "_" $3}'  $outputdir/${n}/$reads_depth_file.MEME._$n.txt > $outputdir/${n}/$reads_depth_file.MEME._$n.list

  #extract each sequences flanking teh modification site
  for seq in $(cat $outputdir/${n}/$reads_depth_file.MEME._$n.list | cut '\t' -f1)
  do
    samtools samtools $ref_genome 
    samtools faidx $ref_genome $seq >> $output_dir/${i}_${n}.fasta
  done

  #rename sequences with both gene ID and position information
  python get_seq.py --ID $outputdir/$reads_depth_file.MEME._$n.list --sequence $outputdir/$reads_depth_file.MEME.seq.fasta

done

#Perform motif enrichemnt
#mkdir $seq
for seq in $(cat ${outputdir}/seq.list | cut -d '.' -f1 )
do
  #mkdir $seq
  cd $seq
  meme ${outputdir}/${seq}.fasta -dna -mod zoops -nmotifs $motif -minw 3 -maxw 15
  mv ${outputdir}/${seq}/meme_out/meme.html ${outputdir}/${seq}.meme.html
done

ELAPSED_TIME_6=$(($SECONDS - $START_TIME_6))
echo "Elapsed time for STEP 6 is" $ELAPSED_TIME_6 "seconds"


START_TIME_7=$SECONDS
echo "******************STEP 6 Performing comparison of coordinates among mods, TEs, gene by sliding windows......******************"

python sliding_window_stats.py --TE $Repeat_gff --MODs ${outputdir}/$reads_depth_file.reads.ratio.csv --Gene $gene_gff --window 10000 --stepsize 2000 -output $window_meta.table

ELAPSED_TIME_7=$(($SECONDS - $START_TIME_7))
echo "Elapsed time for STEP 7 is" $ELAPSED_TIME_7 "seconds"