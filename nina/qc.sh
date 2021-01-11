#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans

if [ $1 == mummer ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/mummer

	nucmer \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/genomes/$i.flye.fasta \
	    $datadir/$i/genomes/$i.canu.fasta
	
	mummerplot \
	    --filter --fat --postscript \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/mummer/$i.delta
	
	mummerplot \
	    --filter --fat --png \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/mummer/$i.delta

	dnadiff \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/genomes/$i.flye.fasta \
	    $datadir/$i/genomes/$i.canu.fasta
    done
fi
