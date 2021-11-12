#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/zymo/truth

if [ $1 == getinfo ] ; then

    awk -F , '{print $1, $12}' $datadir/pacbio/raw/sra_metadata.csv > $datadir/pacbio/raw/pacbio_names.txt

fi

if [ $1 == rename ] ; then
    while read p ; do
	acc=`echo $p | cut -d ' ' -f 1`
	name=`echo $p | cut -d ' ' -f 2`
	if [ $acc != Run ] ; then
	    mv $datadir/pacbio/raw/$acc $datadir/pacbio/raw/$name
	fi
    done < $datadir/pacbio/raw/pacbio_names.txt

fi
