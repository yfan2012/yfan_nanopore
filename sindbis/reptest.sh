#!/bin/bash

##sequencing/basecalling done on the grid

datadir=/dilithium/Data/Nanopore/sindbis
ref=$datadir/refs/sindbis_jane.fasta

if [ $1 == align ] ; then
    for i in reptest_sinv reptest_sinvab ;
    do
	mkdir -p $datadir/$i/align
	minimap2 -ax splice -uf -k14 -t 36 $ref $datadir/$i/fqs/$i.fq | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/$i/align/$i.sorted.bam
	samtools index $datadir/$i/align/$i.sorted.bam

	samtools view -@ 36 -b -F 0x100 $datadir/$i/align/$i.sorted.bam |
	    samtools sort -@ 36 -o $datadir/$i/align/$i.primary.sorted.bam
	samtools index $datadir/$i/align/$i.primary.sorted.bam
	
	mkdir -p $datadir/$i/cov
	bedtools genomecov -d -ibam $datadir/$i/align/$i.sorted.bam > $datadir/$i/cov/$i.cov
	bedtools genomecov -d -ibam $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.cov
    done
fi

