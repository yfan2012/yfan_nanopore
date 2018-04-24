#!/bin/bash

all_refs=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs
readfqs=/atium/Data/Nanopore/cpowgs/170816_BUCC/170816_BUCC.fq.gz
readfas=/atium/Data/Nanopore/cpowgs/170816_BUCC/170816_BUCC.fa
aligndir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_align


seqtk seq -a $readfqs > $readfas

for ref in $all_refs/*.fa ;
do
    reffile=`echo $ref | cut -d '/' -f 8`
    prefix=`echo $reffile | cut -d '.' -f 1`

    minimap2 -t 12 -ax map10k $ref $readfas | samtools view -bS | samtools sort -o $aligndir/$prefix.sorted.bam
    samtools index $aligndir/$prefix.sorted.bam

    samtools faidx $ref
    bedtools genomecov -ibam $aligndir/$prefix.sorted.bam -g $ref.fai -max 1000 > $aligndir/${prefix}_cov.bed
done
