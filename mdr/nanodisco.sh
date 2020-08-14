#!/bin/bash

##recall: must go to ~/software/nanodisco and run `singularity run nd_env` and then go to /home/yfan/Code/.. to run stuff
datadir=/uru/Data/Nanopore/projects/mdr/MDRstool_16
ref=$datadir/metaflye/pcr/MDRstool_16_pcr.assembly.fasta

if [ $1 == preprocess ] ; then
    for i in native pcr ;
    do
	mkdir -p $datadir/nanodisco/$i/preprocess
	nanodisco preprocess \
		  -p 36 \
		  -f $datadir/called_$i/workspace/ \
		  -o $datadir/nanodisco/preprocess \
		  -s $i \
		  -r $ref
    done
fi

	
if [ $1 == coverage ] ; then
    for i in native pcr
    do
	mkdir -p $datadir/nanodisco/$i/coverage/${i}_pcr
	nanodisco coverage \
		  -b $datadir/nanodisco/preprocess/$i.sorted.bam \
		  -o $datadir/nanodisco/preprocess \
		  -r $ref
    done
fi


if [ $1 == current_diff ] ; then
    for i in $samps
    do
	nanodisco chunk_info -r $ref

	mkdir -p $datadir/nanodisco/diffs
	nanodisco difference \
		  -nj 16 \
		  -nc 4 \
		  -p 2 \
		  -w pcr \
		  -n native \
		  -r $datadir/pcr/MDRstool_16_pcr.assembly.fasta \
		  -i $datadir/nanodisco/preprocess \
		  -o $datadir/nanodisco/diffs
		  
    done
fi
	
