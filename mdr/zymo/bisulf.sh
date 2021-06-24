#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/zymo
rawdir=/uru/Data/Nanopore/projects/read_class/zymo
refdir=$rawdir/ref
ref=$refdir/zymo_all.fa

if [ $1 == index ] ; then
    bismark_genome_preparation \
	--parallel 36 \
	$refdir
fi

if [ $1 == bismark ] ; then
    mkdir -p $datadir/bisulfite
    for i in $rawdir/raw/*_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`
	mkdir -p $datadir/bisulfite/$prefix
	bismark \
	    $refdir \
	    -1 $rawdir/raw/${prefix}_1.fastq.gz \
	    -2 $rawdir/raw/${prefix}_2.fastq.gz \
	    -o $datadir/bisulfite/$prefix \
	    -B $prefix \
	    -p 36
    done
fi
    
if [ $1 == bamtosam ] ; then
    for i in $rawdir/raw/*_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`
	samtools view \
		 -@ 36 \
		 -O SAM \
		 -o $datadir/bisulfite/$prefix/${prefix}_pe.sam \
		 $datadir/bisulfite/$prefix/${prefix}_pe.bam
    done
fi

if [ $1 == getmeth ] ; then
    for i in $rawdir/raw/*_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`
	bismark_methylation_extractor \
	    -o $datadir/bisulfite/$prefix \
	    --paired-end \
	    --parallel 18 \
	    --cytosine_report \
	    --CX \
	    --genome_folder $refdir \
	    $datadir/bisulfite/$prefix/${prefix}_pe.sam
	gunzip $datadir/bisulfite/$prefix/*.gz
    done
fi
