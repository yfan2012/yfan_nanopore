#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans


if [ $1 == mummer_mito ] ; then
    mito=$datadir/ref/LProlificans_mito_v1.0.fa
    for i in st5317 ;
    do
	mkdir -p $datadir/$i/final/mummer_mito
	genome=$datadir/$i/final/$i.final.fasta ;
	
	prefix=st5317_correction
	nucmer \
	    -p $datadir/$i/mummer_mito/$prefix \
	    $genome \
	    $mito

	mummerplot \
	    --filter --fat --postscript \
	    -p $datadir/$i/mummer_mito/$prefix \
	    $datadir/$i/mummer_mito/$prefix.delta
	
	mummerplot \
	    --filter --fat --png \
	    -p $datadir/$i/mummer_mito/$prefix \
	    $datadir/$i/mummer_mito/$prefix.delta
	
	dnadiff \
	    -p $datadir/$i/mummer_mito/$prefix \
	    $genome \
	    $mito
    done
fi
