#!/bin/bash

prefix=200708_mdr_stool16native
projdir=/mithril/Data/Nanopore/projects/methbin
clindir=$projdir/mdr
datadir=$projdir/paperfigs/contig_level

barcodes=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes.txt
asm=$clindir/medaka/consensus.fasta

if [ $1 == filter_meth ] ; then
    ##filter out barcode motif related positions in per read meth calls
    python ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $asm \
	   -b $barcodes \
	   -m $clindir/megalodon/${prefix}_asm/per_read_modified_base_calls.txt \
	   -i $clindir/megalodon/${prefix}_asm/per_read_modified_base_calls.txt.idx \
	   -o $datadir/clin_barocdes_methprobs.csv \
	   -t 12
fi
