#!/bin/bash

ssddir=~/data/mdr/disco
datadir=/mithril/Data/Nanopore/projects/methbin/disco

samps='MinION_BF_NAT MinION_CP_NAT MinION_HP_NAT MinION_MH_NAT MinION_NG_NAT MinION_TP_NAT '

ref=$datadir/ref/disco_refs.fasta

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
    for i in $samps ;
    do
	mkdir -p $ssddir/called/$i
	guppy_basecaller \
	    -i $ssddir/raw/$i \
	    -s $ssddir/called/$i \
	    --recursive \
	    --compress_fastq \
	    --flowcell FLO-MIN106 --kit SQK-LSK108 \
	    --device 'cuda:0'

    done
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

if [ $1 == gather ] ; then
    mkdir -p $datadir/fastqs
    for i in $samps ;
    do
	mkdir -p $datadir/fastqs/$i
	cat $datadir/called/$i/pass/*.fastq.gz > $datadir/fastqs/$i/$i.fq.gz
    done
fi


if [ $1 == align ] ; then
    mkdir -p $datadir/align
    for i in $samps ;
    do
	mkdir -p $datadir/align/$i
	fq=$datadir/fastqs/$i/$i.fq.gz
	minimap2 -t 36 -x map-ont $ref $fq \
		 > $datadir/align/$i/$i.paf
    done
fi



if [ $1 == filter ] ; then
    for i in $samps ;
    do
	for j in 1 5 10 ;
	do
	    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
		   -a $datadir/align/$i/$i.paf \
		   -r $ref \
		   -o $datadir/barcode/$i/${i}_barcodes_filtered.$j.txt \
		   -m $datadir/barcode/$i/${i}_barcodes.txt \
		   -b ~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt \
		   -q 40 \
		   -l 5000 \
		   -n $j \
		   -v
	done
    done
fi




















