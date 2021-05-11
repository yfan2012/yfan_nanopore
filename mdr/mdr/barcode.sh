#!/bin/bash

ssddir=~/data/mdr/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/disco
prefix=200708_mdr_stool16native

if [ $1 == megaidx ] ; then
    for i in $samps ;
    do
	python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
		-i $ssddir/megalodon/$i/per_read_modified_base_calls.txt \
		-o $ssddir/megalodon/$i/per_read_modified_base_calls.txt.idx
    done
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

if [ $1 == barcode ] ; then
    mkdir -p $ssddir/barcode
    for i in $samps ;
    do
	mkdir -p $ssddir/barcode/$i
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $ssddir/megalodon/$i/per_read_modified_base_calls.txt \
               -i $ssddir/megalodon/$i/per_read_modified_base_calls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt \
               -o $ssddir/barcode/$i/${i}_barcodes.txt \
	       -n 100000000000 \
               -t 36 ;} &> $ssddir/barcode/$i/${i}_time.txt
    done
fi


if [ $1 == copy ] ; then
    cp -r $ssddir/barcode $datadir/
    cp -r $ssddir/called $datadir/
    cp -r $ssddir/megalodon $datadir/
fi
