#!/bin/bash

prefix=20190809_zymo_control
datadir=/mithril/Data/Nanopore/projects/methbin/zymo
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == index ] ; then
    ##second indexing
    ##moved old index to per_read_modified_base_calls.txt.idx_old
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
	    -o $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx
fi


if [ $1 == call_test ] ; then
    mkdir -p $datadir/barcode
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.small.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
           -o $datadir/barcode/test.txt \
           -t 12 ;} &> $datadir/barcode/test_time.txt
fi

	    
if [ $1 == call ] ; then
    mkdir -p $datadir/barcode_v2
    mkdir -p $datadir/barcode_v2/$prefix
    
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
           -o $datadir/barcode_v2/$prefix/${prefix}_barcodes20.txt \
	   -n 100000000000 \
           -t 36 ;} &> $datadir/barcode_v2/$prefix/${prefix}_time20.txt
fi


if [ $1 == filtcommon ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -a $datadir/align/$prefix.paf \
	   -r $ref \
	   -o $datadir/barcode_v2/$prefix/${prefix}_motifcounts20.txt \
	   -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
	   -q 40

fi
