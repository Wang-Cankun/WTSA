#!/bin/bash
echo "1.bowtie2 mapping"
module add gcc
module add bowtie2
module add samtools
bowtie2 -p 8 -x /pylon5/ci4s87p/cankun/genome/Esche.fasta -f /pylon5/ci4s87p/cankun/fasta/$1.fasta -S /pylon5/ci4s87p/cankun/fasta/$1.sam
echo "bam sorting"
samtools view -@ 8 -Sb  /pylon5/ci4s87p/cankun/fasta/$1.sam  >  /pylon5/ci4s87p/cankun/fasta/$1.bam
samtools sort -@ 8 -m 10G /pylon5/ci4s87p/cankun/fasta/$1.bam -o /pylon5/ci4s87p/cankun/fasta/$1.sorted.bam
echo "bam indexing"
samtools index /pylon5/ci4s87p/cankun/fasta/$1.sorted.bam 
echo "bedgraph creating"
#bedtools genomecov -ibam /pylon5/ci4s87p/cankun/fasta/$1.sorted.bam -bg > $1.bg

