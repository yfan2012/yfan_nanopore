#!/bin/bash

##code to look at dunlop pilon polishing
datadir=/kyber/Data/Nanopore/projects/asmpolish/dunlop/180714_dunlop

if [ $1 == assemble ] ; then
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	mkdir -p $datadir/$i/${i}_assemble
	canu \
	    -p $i -d $datadir/$i/${i}_assemble \
	    genomeSize=4.8m \
	    -nanopore-raw $datadir/$i/$i.fastq
    done
fi

   

