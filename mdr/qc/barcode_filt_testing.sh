#!/bin/bash

#datadir=~/data/mdr/qc
datadir=/mithril/Data/Nanopore/projects/methbin
ref=/mithril/Data/Nanopore/projects/methbin/reference/allsamps.fa


if [ $1 == test_thresholds ] ; then
    ##try setting really stricter thresholds to see if barcodes resolve better
    ##know that cmod roc sets thresh around 1 and amod sets thresh at around 0
    #for i in neb15 neb17 neb19 nebdcm neb11 ;

    ##for i in nebdcm ;
    for i in neb15 neb17 neb19 neb11 ;
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

if [ $1 == test ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode.py \
	   -m $datadir/megalodon/test/testdcm_mod_basecalls.txt \
	   -i $datadir/megalodon/test/testdcm_mod_basecalls.txt.idx \
	   -r $ref \
	   -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
	   -o $datadir/megalodon/test/testdcm_barcodes.txt \
	   -t 12
fi
	    
	    
if [ $1 == test_thresh ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode.py \
	   -m $datadir/megalodon/test/testdcm_mod_basecalls.txt \
	   -i $datadir/megalodon/test/testdcm_mod_basecalls.txt.idx \
	   -r $ref \
	   -a 0 \
	   -c 1 \
	   -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
	   -o $datadir/megalodon/test/testdcm_barcodes.txt \
	   -t 12
fi
