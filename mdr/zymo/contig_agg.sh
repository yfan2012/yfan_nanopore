#!/bin/bash

projdir=/mithril/Data/Nanopore/projects/methbin/zymo
datadir=$projdir/contig_agg
prefix=20190809_zymo_control
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == test ] ; then
    mkdir -p $datadir
    mkdir -p $datadir/test

    time ( python ~/Code/yfan_meth/utils/megalodon_agg_experimental.py \
	   -m $projdir/megalodon/$prefix/per_read_modified_base_calls.txt \
	   -i $projdir/megalodon/$prefix/per_read_modified_base_calls.txt.small.idx \
	   -r $ref \
	   -o $datadir/test/test.txt \
	   -t 12 \
	   -v )
fi

if [ $1 == agg ] ; then
    mkdir -p $datadir
    mkdir -p $datadir/$prefix

    python ~/Code/yfan_meth/utils/megalodon_agg.py \
	   -m $projdir/megalodon/$prefix/per_read_modified_base_calls.txt \
	   -i $projdir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
	   -r $ref \
	   -o $datadir/$prefix/$prefix.meth_report.txt \
	   -t 12 \
	   -v
fi

