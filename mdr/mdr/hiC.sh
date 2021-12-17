#!/bin/bash

projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/mdr/hiC
prefix=200708_mdr_stool16native

catdir=/atium/Data/ref/CATdb
if [ $1 == BAT ] ; then
    mkdir -p $datadir/bin_id
    mkdir -p $datadir/bin_id/BAT
    
    CAT bins \
	-b $datadir/clusters \
	-d $catdir/db \
	-t $catdir/tax \
	-o $datadir/bin_id/BAT/$prefix.BAT \
	-s fasta \
	--force \
	--sensitive
fi


if [ $1 == namebins ] ; then
    CAT add_names \
        -i $datadir/bin_id/BAT/$prefix.BAT.bin2classification.txt \
        -t $catdir/tax \
        -o $datadir/bin_id/BAT/$prefix.BAT.names.txt
fi
