#!/bin/bash

datadir=~/data/mdr/qc
ref=/mithril/Data/Nanopore/projects/methbin/reference/allsamps.fa

if [ $1 == test ] ; then
    mkdir -p $datadir/barcode
    for i in test ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
	       -m $datadir/test.txt \
	       -i $datadir/test.txt.idx \
	       -r $ref \
	       -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
	       -o $datadir/barcode/test.txt \
	       -t 12 ;} &> $datadir/barcode/test_time.txt
    done
fi


    
if [ $1 == call ] ; then
    for i in neb15 neb17 neb19 nebdcm ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt \
               -i $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
               -o $datadir/barcode/${i}_barcodes.txt \
               -t 36 ;} &> $datadir/barcode/${i}_time.txt
    done
fi

	
