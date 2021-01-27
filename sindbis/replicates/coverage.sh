#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis/replicates
ref=/dilithium/Data/Nanopore/sindbis/refs/sindbis_jane.fasta
i=$2

if [ $1 == align ] ; then

    mkdir -p $datadir/$i/align
    minimap2 -a -uf -k14 -t 36 $ref $datadir/$i/fqs/$i.fq.gz | \
        samtools view -@ 36 -b | \
        samtools sort -@ 36 -o $datadir/$i/align/$i.sorted.bam
    samtools index $datadir/$i/align/$i.sorted.bam
    
    samtools view -@ 36 -b -F 0x100 $datadir/$i/align/$i.sorted.bam |
        samtools sort -@ 36 -o $datadir/$i/align/$i.primary.sorted.bam
    samtools index $datadir/$i/align/$i.primary.sorted.bam
    
fi


if [ $1 == cov ] ; then
    mkdir -p $datadir/$i/cov
    bedtools coverage -d -a $ref.bed -b $datadir/$i/align/$i.sorted.bam > $datadir/$i/cov/$i.cov
    bedtools coverage -d -a $ref.bed -b $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.cov
fi

if [ $1 == genomecov ] ; then
    mkdir -p $datadir/$i/cov
    bedtools genomecov -d -ibam $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.genomecov
fi