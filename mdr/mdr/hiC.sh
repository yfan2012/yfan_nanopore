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

if [ $1 == BAT_single ] ; then
    ##downstream naming doesn't tolerate two possible classifications
    ##force single classification using the f parameter
    mkdir -p $datadir/bin_id
    mkdir -p $datadir/bin_id/BAT_single
    
    CAT bins \
	-b $datadir/clusters \
	-d $catdir/db \
	-t $catdir/tax \
	-o $datadir/bin_id/BAT_single/$prefix.BAT \
	-s fasta \
	-f .5 \
	--force \
	--sensitive
fi


if [ $1 == namebins ] ; then
    CAT add_names \
        -i $datadir/bin_id/BAT/$prefix.BAT.bin2classification.txt \
        -t $catdir/tax \
        -o $datadir/bin_id/BAT/$prefix.BAT.names.txt
fi


if [ $1 == summarise ] ; then
    CAT add_names \
	-i $datadir/bin_id/BAT_single/$prefix.BAT.bin2classification.txt \
	-t $catdir/tax \
	-o $datadir/bin_id/BAT_single/$prefix.BAT.names_official.txt \
	--only_official
    CAT summarise \
	-i $datadir/bin_id/BAT_single/$prefix.BAT.names_official.txt \
	-o $datadir/bin_id/BAT_single/$prefix.BAT.summary.txt
fi
    

if [ $1 == abricate ] ; then
    mkdir -p $projdir/mdr/amr

    binfile=$datadir/$prefix.hiC.fasta

    for i in ecoh card ncbi resfinder plasmidfinder vfdb ecoli_vf megares argannot ;
    do
	abricate \
	    --threads 36 \
	    --db $i \
	    $binfile > $projdir/mdr/amr/$prefix.hiC.$i.tsv
    done
fi
    
    
