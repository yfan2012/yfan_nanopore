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
    
