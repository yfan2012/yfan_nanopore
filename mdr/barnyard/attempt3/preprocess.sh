#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/allrefs.fa

prefix=211005_mdr_barnyard_mix

if [ $1 == basecall ] ; then
    mkdir -p $ssddir/called

    for i in 5 6 ;
    do
	samp=${prefix}${i}
	
	mkdir -p $ssddir/called/$samp
	guppy_basecaller \
	    -i $ssddir/raw/$samp/no_sample/*/fast5_pass \
	    -s $ssddir/called/$samp \
	    --recursive \
	    --compress_fastq \
	    --flowcell FLO-FLG001 --kit SQK-RAD004 \
	    --device 'cuda:0'
    done
fi

if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs

    for i in 5 6 ;
    do
	samp=${prefix}${i}
	cat $ssddir/called/$samp/pass/*fastq.gz > $datadir/fastqs/$samp.fastq.gz
    done
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/align
    for i in 5 6 ;
    do
	samp=${prefix}${i}
	
	mkdir -p $datadir/align/$samp
	fq=$datadir/fastqs/$samp.fastq.gz

	minimap2 -t 36 -x map-ont $ref $fq \
		 > $datadir/align/$samp.paf
    done
fi

if [ $1 == alignbam ] ; then

    for i in 5 6 ;
    do
	samp=${prefix}${i}
	fq=$datadir/fastqs/$samp.fastq.gz

	minimap2 -t 36 -ax map-ont $ref $fq | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/align/$samp.sorted.bam
	samtools index $datadir/align/$samp.sorted.bam
		 
    done
fi

	
