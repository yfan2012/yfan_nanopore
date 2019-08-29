#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/mdr5

if [ $1 == transfer ] ; then
    for i in $datadir/* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	echo $prefix
	scp $i/$prefix.fastq whack:/uru/Data/Nanopore/projects/mdr_phase/data/$prefix/
    done
fi

if [ $1 == gather ] ; then
    for i in $datadir/* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	cat $i/called/*/workspace/pass/*fastq > $i/$prefix.fastq
    done
fi
