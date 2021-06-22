#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/

if [ $1 == bisulf ] ; then
    for i in SRR8137538 SRR8137545 SRR8137540 SRR8137541 SRR8137542 SRR8137547 SRR8137544 ;
    do
	fastq-dump \
	    -O $rawdir/zymo/raw/ \
	    --gzip \
	    --split-files \
	    $i &
    done
fi

if [ $1 == rename ] ; then
    rename 's/SRR8137538/SAEN/' $rawdir/zymo/raw/*.fastq.gz
    rename 's/SRR8137541/PSAE/' $rawdir/zymo/raw/*.fastq.gz
    rename 's/SRR8137542/LIMO/' $rawdir/zymo/raw/*.fastq.gz
    rename 's/SRR8137544/ENFA/' $rawdir/zymo/raw/*.fastq.gz
    rename 's/SRR8137545/ESCO/' $rawdir/zymo/raw/*.fastq.gz
    rename 's/SRR8137547/BASU/' $rawdir/zymo/raw/*.fastq.gz
    rename 's/SRR8137540/STAU/' $rawdir/zymo/raw/*.fastq.gz
fi
