#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/allrefs.fa

##prefix=210908_mdr_barnyard_mix1
##prefix=210908_mdr_barnyard_mix2
##prefix=210912_mdr_barnyard_mix3
prefix=210912_mdr_barnyard_mix4

if [ $1 == basecall ] ; then
    mkdir -p $ssddir/called
    mkdir -p $ssddir/called/$prefix
    
    guppy_basecaller \
	-i $ssddir/raw/$prefix/no_sample/*/fast5_pass \
	-s $ssddir/called/$prefix \
	--recursive \
	--compress_fastq \
	--flowcell FLO-FLG001 --kit SQK-RAD004 \
	--device 'cuda:0'
fi

if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs

    cat $ssddir/called/$prefix/pass/*fastq.gz > $datadir/fastqs/$prefix.fastq.gz
fi

##name=210908_mdr_barnyard_mix
name=210912_mdr_barnyard_mix
if [ $1 == align ] ; then
    mkdir -p $datadir/align

    ##for i in 1 2 ;
    for i in 3 4 ;
    do
	prefix=${name}$i
	fq=$datadir/fastqs/$prefix.fastq.gz

	minimap2 -t 36 -x map-ont $ref $fq \
		 > $datadir/align/$prefix.paf
    done
fi

if [ $1 == alignbam ] ; then

    ##for i in 1 2 ;
    for i in 3 4 ;
    do
	prefix=${name}$i
	fq=$datadir/fastqs/$prefix.fastq.gz

	minimap2 -t 36 -ax map-ont $ref $fq | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/align/$prefix.sorted.bam
	samtools index $datadir/align/$prefix.sorted.bam
		 
    done
fi

	
