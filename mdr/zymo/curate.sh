#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/zymo/raw
ssddir=~/data/mdr/zymo
prefix=20190809_zymo_control
datadir=/mithril/Data/Nanopore/projects/methbin


ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == call ] ; then
    mkdir -p $datadir/zymo/barcode
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/zymo/barcodes_zymo_curated.txt \
           -o $datadir/zymo/barcode/${prefix}_curated.txt \
	   -n 100000000000 \
           -t 12 ;} &> $datadir/zymo/barcode/${prefix}_curated.time.txt
fi

if [ $1 == filtcommon ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -a $datadir/zymo/align/$prefix.paf \
	   -r $ref \
	   -o $datadir/zymo/barcode/${prefix}_motifcounts_curated.txt \
	   -m $datadir/zymo/barcode/${prefix}_barcodes_curated.txt \
	   -b ~/Code/yfan_nanopore/mdr/zymo/barcodes_zymo_curated.txt \
	   -q 40 \
	   -v
fi


