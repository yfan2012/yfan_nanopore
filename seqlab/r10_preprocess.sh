#!/bin/bash

##assemble on whack
prefix=190706_nivar_r10
datadir=/kyber/Data/Nanopore/Analysis/$prefix

##called using guppy 3.2.1 in docker container


if [ $1 == fqs ] ; then
    cat $datadir/called/*fastq > $datadir/$prefix.fq
fi

    
if [ $1 == longfqs ] ; then
    python2 ~/Code/utils/fastq_long.py -i $datadir/$prefix.fq -o $datadir/${prefix}_3k.fq -l 3000
fi


if [ $1 == assemble ] ; then
    ##kyber might get disconnected tonight, so I copied fqs to whack on board storage to run canu
    datadir=~/data/nivar
    canu \
	-p nivar -d $datadir/canu \
	genomeSize=15m \
	-nanopore-raw $datadir/${prefix}_3k.fq
fi

if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed
    for i in cDNA gDNA ;
    do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 72 -phred33 \
	     $datadir/${i}_illumina/CANI_${i}_R1.fastq.gz $datadir/${i}_illumina/CANI_${i}_R2.fastq.gz \
	     $datadir/trimmed/CANI_${i}_forward_paired.fq.gz $datadir/trimmed/CANI_${i}_forward_unpaired.fq.gz \
	     $datadir/trimmed/CANI_${i}_reverse_paired.fq.gz $datadir/trimmed/CANI_${i}_reverse_unpaired.fq.gz \
	     ILLUMINACLIP:idt_ud_indexes.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi
