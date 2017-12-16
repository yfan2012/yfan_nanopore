#!/bin/bash

all_refs=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs
reads=/atium/Data/NGS/Raw/171002_bucc
aligndir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_align_ill

for ref in $all_refs/*.fa ;
do
    prefix=${ref%.fa}
    refpre=`echo $prefix | cut -d '/' -f 8`
    bowtie2-build $ref $prefix
    bowtie2 -p 12 -x $prefix -1 $reads/bucc_S1_L001_R1_001.fastq.gz -2 $reads/bucc_S1_L001_R2_001.fastq.gz | samtools view -bS | samtools sort -o $aligndir/$refpre.sorted.bam
    samtools index $aligndir/$refpre.sorted.bam

    ##index already made from nanopore alignment
    ##samtools faidx $ref
    bedtools genomecov -ibam $aligndir/$refpre.sorted.bam -g $ref.fai -max 1000 > $aligndir/${refpre}_cov.bed
    
done
