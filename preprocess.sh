#/bin/bash

PEAK_FILE=$1
BAM_FILE=$2
REFERENCE_FILE=$3
OUTPUTFILE=$4

bedtools coverage -a $PEAK_FILE -b $BAM_FILE -d | cut -f4,12 | awk -f group.awk > tmp.coverage.txt

bedtools getfasta -fi $REFERENCE_FILE -bed $PEAK_FILE  > tmp.fasta

awk '(NR)%2==0 {getline this<"tmp.coverage.txt";print this} 1' tmp.fasta | tac | sed '/N/,+2d' | tac > $OUTPUTFILE

#rm tmp.coverage.txt
#rm tmp.fasta

