#!/bin/bash

##six reads that span the genome according to duncan
##data is in the email
datadir=/uru/Data/Nanopore/projects/circ

if [ $1 == ecoli6_canu ] ; then
    ##try canu
    
    canu \
	-p ecoli6 \
	-d $datadir/ecoli6/canu \
	genomeSize=4.5m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/ecoli6/ecoli_6reads.fq
fi

if [ $1 == ecoli6_flye ] ; then
    ##try flye

    flye \
	--nano-raw $datadir/ecoli6/ecoli_6reads.fq \
	-g 1m \
	-o $datadir/ecoli6/flye \
	-t 24
fi
    
       
   
    
