#!/bin/bash

ssddir=~/data/mdr/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/mdr
prefix=200708_mdr_stool16native

ref=$datadir/ref/mdr_refs.fa

if [ $1 == megaidx ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt \
	    -o $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt.idx
fi


if [ $1 == guppy_call ] ; then
    mkdir -p $ssddir/called

    guppy_basecaller \
	-i $ssddir/raw \
	-s $ssddir/called \
	--recursive \
	--compress_fastq \
	--flowcell FLO-MIN106 --kit SQK-LSK109 \
	--device 'cuda:0'
    
fi

if [ $1 == copy ] ; then
    cp -r $ssddir/called $datadir/
    cp -r $ssddir/megalodon $datadir/
fi


if [ $1 == barcode ] ; then
    mkdir -p $datadir/barcode
    mkdir -p $datadir/barcode/$prefix
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
           -o $datadir/barcode/$prefix/${prefix}_barcodes.txt \
           -t 12 ;} &> $datadir/barcode/$prefix/${prefix}_time.txt
fi


