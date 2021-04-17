#!/bin/bash

#datadir=~/data/mdr/qc
datadir=/mithril/Data/Nanopore/projects/methbin
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
    for i in neb15 neb17 neb19 nebdcm neb11 ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt \
               -i $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
               -o $datadir/barcode/${i}_barcodes.txt \
               -t 72 ;} &> $datadir/barcode/${i}_time.txt
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
               -o $datadir/barcode/${i}_barcodes.txt \
               -t 36 ;} &> $datadir/barcode/${i}_time.txt
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
    ##this didn't really work that well. forget this.
fi


if [ $1 == test_thresholds ] ; then
    ##try setting really stricter thresholds to see if barcodes resolve better
    ##know that cmod roc sets thresh around 1 and amod sets thresh at around 0
    #for i in neb15 neb17 neb19 nebdcm neb11 ;

    for i in nebdcm ;
    do
	for j in 9 8 7 6 5 4 3 2 1 ;
	do
	    cbound=1.$j
	    abound=.$j
	    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
		   -m $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt \
		   -i $datadir/megalodon/$i/$i/${i}_mod_basecalls.txt.idx \
		   -r $ref \
		   -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
		   -a $abound \
		   -c $cbound \
		   -o $datadir/barcode/${i}_barcodes_${j}.txt \
		   -t 36 ;} &> $datadir/barcode/${i}_time_${j}.txt
	    
	done
    done
fi

	    


	    
