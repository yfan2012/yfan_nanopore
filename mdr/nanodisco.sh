#!/bin/bash

##recall: must go to ~/software/nanodisco and run `singularity run nd_env` and then go to /home/yfan/Code/.. to run stuff
datadir=/uru/Data/Nanopore/projects/mdr/MDRstool_16


if [ $1 == preprocess ] ; then
    for i in native pcr ;
    do
	mkdir -p $datadir/nanodisco/$i/preprocess
	nanodisco preprocess \
		  -p 36 \
		  -f $datadir/called_$i/workspace/ \
		  -s $i \
		  -o $datadir/$i/nanodisco/preprocess \
		  -r $datadir/metaflye/${i}_100m/MDRstool_16_${i}_100m.assembly.fasta
    done
fi

	
if [ $1 == coverage ] ; then
    for i in $samps
    do
	mkdir -p $datadir/$i/nanodisco/coverage
	nanodisco coverage \
		  -b $datadir/$i/nanodisco/preprocess/$i.sorted.bam \
		  -r $datadir/$i/metaflye/150m/assembly.fasta \
		  -o $datadir/$i/nanodisco/coverage
    done
fi


if [ $1 == current_diff ] ; then
    for i in $samps
    do
	nanodisco chunk_info -r $datadir/$i/metaflye/150m/assembly.fasta

	nanodisco difference \
		  -nj 16 \
		  -nc 4 \
		  -p 2 \
		  -r $datadir/$i/metaflye/150m/assembly.fasta \
		  -i $datadir/$i/nanodisco/preprocess 
		  
    done
fi
	
