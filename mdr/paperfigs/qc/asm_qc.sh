#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin
mdrdir=$datadir/mdr

dbxdir=~/gdrive/mdr/paperfigs/qc

if [ $1 == hic ] ; then
    python3 ~/Code/utils/qc/asm_assess.py \
	    -i $mdrdir/hiC/200708_mdr_stool16native.hiC.fasta \
	    -p mdrhic >> ~/gdrive/mdr/paperfigs/figs/asmstats.csv
    
fi
