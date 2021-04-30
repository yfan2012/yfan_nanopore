#!/bin/bash

##getting data organized for qc of megalodon/guppy meth calling.
##mostly stolen from yfan_meth/rerio
##subset stuff is already set up pretty nicely in mithril

rawdir=/mithril/Data/Nanopore/projects/methbin
datadir=~/data/mdr/qc

if [ $1 == cp_multiraw_sub ] ; then
    mkdir -p $datadir
    cp -r $rawdir/multiraw_sub $datadir
fi

if [ $1 == clean_ctrl ] ; then
    for i in neb12 neb13 neb15 neb16 neb17 neb19 nebdcm ;
    do
	rm $datadir/multiraw_sub/$i/neb11*
    done
fi


if [ $1 == basecall ] ; then
    ##need basecalls and alignment to filter barcoded reads based on megalodon calls
    mkdir -p $datadir/fastqs
    for i in neb11 neb12 neb13 neb15 neb16 neb17 neb19 nebdcm ;
    do
	guppy_basecaller \
	    -i $datadir/multiraw_sub/$i \
	    -s $datadir/fastqs/$i \
	    --compress_fastq \
	    --flowcell FLO-MIN106 --kit SQK-LSK108 \
	    --device 'cuda:0'
    done
fi


if [ $1 == move_fastqs ] ; then
    for i in neb11 neb12 neb13 neb15 neb16 neb17 neb19 nebdcm ;
    do
	cat $datadir/fastqs/$i/pass/*fastq.gz > $rawdir/fastqs/$i/$i.fq.gz
    done
fi

ref=$rawdir/reference/allsamps.fa
if [ $1 == align ]  ; then
    for i in neb11 neb12 neb13 neb15 neb16 neb17 neb19 nebdcm ;
    do
	fq=$rawdir/fastqs/$i/$i.fq.gz
	minimap2 -t 36 -x map-ont $ref $fq \
		 > $rawdir/align/$i/${i}_sub.paf
    done
fi
