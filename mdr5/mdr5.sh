#!/bin/bash

datadir=/uru/Data/Nanopore/projects/mdr_phase/data

if [ $1 == meta_asm ] ; then
    for i in recswab_1r mdrstool_16 mdrstool_19 ;
    do
	mkdir -p $datadir/$i/flye
	flye --nano-raw $datadir/$i/$i.fastq --meta -g 150m -o $datadir/$i/flye -i 5 -t 54
    done
fi
