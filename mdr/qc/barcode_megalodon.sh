#!/bin/bash

#datadir=~/data/mdr/qc/megalodon
datadir=/mithril/Data/Nanopore/projects/methbin
samps='neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm'

if [ $1 == index_megalodon_output ] ; then
    for i in $samps ;
    do
	time (python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
		-i $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt \
		-o $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt.idx ) \
	    &> $datadir/megalodon/$i/$i/${i}_mod_index_time.txt
    done
fi



ref=/mithril/Data/Nanopore/projects/methbin/reference/allsamps.fa

if [ $1 == test ] ; then
    for i in test ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
	       -m $datadir/megalodon/test/testdcm_mod_basecalls.txt \
	       -i $datadir/megalodon/test/testdcm_mod_basecalls.txt.idx \
	       -r $ref \
	       -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
	       -o $datadir/megalodon/test/testdcm_barcodes.txt \
	       -t 12 ;} &> $datadir/megalodon/test/testdcm_time.txt
    done
fi


    
if [ $1 == call ] ; then
    for i in neb15 neb17 neb19 nebdcm neb11 ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt \
               -i $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
               -o $datadir/barcode/${i}_barcodes.txt \
               -t 72 ;} &> $datadir/barcode/qc/${i}_time.txt
    done
fi

if [ $1 == call_unmeth ] ; then
    for i in neb11 ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt \
               -i $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
               -o $datadir/barcode/qc/${i}_barcodes.txt \
               -t 36 ;} &> $datadir/barcode/qc/${i}_time.txt
    done
fi



if [ $1 == test_agg ] ; then
    mkdir -p $datadir/barcode_agg
    for i in test ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_aggregate_barcode.py \
	       -m $datadir/test.txt \
	       -i $datadir/test.txt.idx \
	       -r $ref \
	       -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
	       -o $datadir/barcode_agg/test.txt \
	       -n 10000 \
	       -t 12 ;} &> $datadir/barcode_agg/test_time.txt
    done
fi


if [ $1 == call_agg ] ; then
    ##for i in neb15 neb17 neb19 nebdcm ;
    for i in neb19 nebdcm ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_aggregate_barcode.py \
               -m $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt \
               -i $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
               -o $datadir/barcode_agg/${i}_barcodes.txt \
	       -n 5000 \
               -t 72 ;} &> $datadir/barcode_agg/${i}_time.txt
    done
    ##this didn't really work that well. forget it.
fi


	    
