#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans/rna

if [ $1 == concat ] ; then
    mkdir -p $datadir/reads

    cat $datadir/raw/*1.fq.gz > $datadir/reads/5317_R1.fq.gz
    cat $datadir/raw/*2.fq.gz > $datadir/reads/5317_R2.fq.gz
fi

if [ $1 == trim ] ; then
    ##make adapters file. single adapter files don't have trailing newlines...
    cat ~/software/Trimmomatic-0.39/adapters/*fa > adapters.fa
    ##put in appropriate newlines
    sed -i -e 's/>/\'$'\n>/g' adapters.fa
    ##get rid of empty lines
    sed -i '/^$/d' adapters.fa 
    
    mkdir -p $datadir/trimmed
    
    java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	 -threads 36 -phred33 \
	 $datadir/reads/5317_R1.fq.gz $datadir/$i/reads/5317_R2.fq.gz \
	 $datadir/trimmed/5317_fwd_paired.fq.gz $datadir/trimmed/5317_fwd_unpaired.fq.gz \
	 $datadir/trimmed/5317_rev_paired.fq.gz $datadir/trimmed/5317_rev_unpaired.fq.gz \
	 ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
fi

if [ $1 == trinity ] ; then
    mkdir -p $datadir/trinity

    Trinity \
	--seqType fq \
	--max_memory 250G \
	--CPU 54 \
	--left $datadir/trimmed/5317_fwd_paired.fq.gz \
	--right $datadir/trimmed/5317_rev_paired.fq.gz \
	--output $datadir/trinity
fi
