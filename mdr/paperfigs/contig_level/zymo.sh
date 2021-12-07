#!/bin/bash

prefix=20190809_zymo_control
projdir=/mithril/Data/Nanopore/projects/methbin
zymodir=$projdir/zymo
datadir=$projdir/paperfigs/contig_level

barcodes=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/zymo_barcodes.txt
asm=$zymodir/medaka/consensus.fasta

if [ $1 == filter_meth ] ; then
    ##filter out barcode motif related positions in per read meth calls
    python ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $asm \
	   -b $barcodes \
	   -m $zymodir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt \
	   -i $zymodir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt.idx \
	   -o $datadir/zymo_barocdes_methprobs.csv \
	   -t 12
fi

