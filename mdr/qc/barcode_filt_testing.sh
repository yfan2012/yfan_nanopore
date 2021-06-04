#!/bin/bash

#datadir=~/data/mdr/qc
projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/barcode
ref=/mithril/Data/Nanopore/projects/methbin/reference/allsamps.fa


if [ $1 == test ] ; then
    testdir=$projdir/megalodon/test
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
	   -m $testdir/testdcm_mod_basecalls.txt \
	   -i $testdir/testdcm_mod_basecalls.txt.idx \
	   -r $ref \
	   -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
	   -o $testdir/testdcm_barcodes.txt \
	   -t 12 ;} &> $testdir/testdcm_time.txt
    { time python ~/Code/yfan_meth/utils/megalodon_barcode_test.py \
	   -m $testdir/testdcm_mod_basecalls.txt \
	   -i $testdir/testdcm_mod_basecalls.txt.idx \
	   -r $ref \
	   -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
	   -o $testdir/testdcm_barcodes2.txt \
	   -t 12 ;} &> $testdir/testdcm_time2.txt

fi

if [ $1 == test_threshold_option ] ; then
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

if [ $1 == try_thresholds ] ; then
    ##try setting stricter thresholds to see if barcodes resolve better
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
		   -o $datadir/qc/${i}_barcodes_${j}.txt \
		   -t 36 ;} &> $datadir/qc/${i}_time_${j}.txt
	    
	done
    done
fi

if [ $1 == test_filterer ] ; then
    ##test the code that filters based on alignment info (readlen, mapq, nummotifs, etc)    
    for i in neb15 neb17 neb19 neb11 nebdcm ;
    do
	for j in 5 10 15 20 25 30 ;
	do
	    echo $i
	    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
		   -m $datadir/qc/${i}_barcodes.txt \
		   -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
		   -a $projdir/align/$i/${i}_sub.paf \
		   -r $ref \
		   -q 30 \
		   -l 5000 \
		   -n $j \
		   -o $datadir/qc/${i}_barcodes_filtered_${j}_motifs.txt
	done
    done
fi


	    


