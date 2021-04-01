#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans

if [ $1 == trf ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/repeats
	asm=$datadir/$i/final/$i.final.fasta
	
	~/software/trf $asm 2 7 7 80 10 50 600 -ngs > $datadir/$i/repeats/$i.trf.txt
    done
fi
	
