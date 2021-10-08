#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/allrefs.fa

name=211005_mdr_barnyard_mix

if [ $1 == megaidx ] ; then
    for i in 5 6 ;
    do
	(prefix=${name}$i
	python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
		-i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
		-o $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx ) &
    done
fi
 
if [ $1 == barcode ] ; then
    mkdir -p $datadir/barcode
    for i in 5 6 ;
    do
	prefix=${name}$i
	python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
           -o $datadir/barcode/${prefix}_barcodes.txt \
           -t 12
    done
fi


if [ $1 == motifcounts ] ; then
    for i in 5 6 ;
    do
	prefix=${name}$i
	python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	       -r $ref \
	       -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
	       -a $datadir/align/$prefix.paf \
	       -m $datadir/barcode/${prefix}_barcodes.txt \
	       -o $datadir/barcode/${prefix}_barcodes_motifcounts.txt \
	       -q 40
    done
fi
