#!/bin/bash

prefix=20190809_zymo_control
datadir=/mithril/Data/Nanopore/projects/methbin/zymo
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == test_kmer ] ; then
    mkdir -p $datadir/kmer_control

    { time python ~/Code/yfan_meth/utils/megalodon_barcode_kmer_control.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.small.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
           -o $datadir/kmer_control/test.txt \
           -t 12 ;} &> $datadir/kmer_control/test_time.txt
fi
    
if [ $1 == run_control ] ; then
    mkdir -p $datadir/kmer_control/$prefix

    { time python ~/Code/yfan_meth/utils/megalodon_barcode_kmer_control.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
           -o $datadir/kmer_control/$prefix/${prefix}_barcodes20.txt \
           -t 12 ;} &> $datadir/kmer_control/$prefix/${prefix}_time20.txt
fi


if [ $1 == test_motifs ] ; then
    mkdir -p $datadir/kmer_control

    { time python ~/Code/yfan_meth/utils/megalodon_barcode_motif_control.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.small.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
           -o $datadir/kmer_control/test.txt \
           -t 12 ;} &> $datadir/kmer_control/test_time.txt
fi

if [ $1 == run_motifs ] ; then
    mkdir -p $datadir/kmer_control/${prefix}_motif

    { time python ~/Code/yfan_meth/utils/megalodon_barcode_motif_control.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt \
           -o $datadir/kmer_control/${prefix}_motif/${prefix}_barcodes20.txt \
           -t 12 ;} &> $datadir/kmer_control/${prefix}_motif/${prefix}_time20.txt
fi
