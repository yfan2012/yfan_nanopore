#!/bin/bash

ref=/atium/Data/Nanopore/cpowgs/bucc/170816_BUCC.contigs.fasta
bucc=/atium/Data/Nanopore/cpowgs/bucc
readfqs=$bucc/170816_BUCC.fq.gz
readfas=$bucc/170816_BUCC.fa

seqtk seq -a $readfqs > $readfas

minimap2 -t 12 -ax map10k $ref $readfas > $bucc/170816_BUCC.sam

samtools view -bS $bucc/170816_BUCC.sam > $bucc/170816_BUCC.bam
samtools sort -o $bucc/170816_BUCC.sorted.bam $bucc/170816_BUCC.bam
samtools index $bucc/170816_BUCC.sorted.bam

samtools faidx $ref
bedtools genomecov -ibam $bucc/170816_BUCC.sorted.bam -g $bucc/170816_BUCC.contigs.fasta.fai -max 1000 > $bucc/cov.bed
