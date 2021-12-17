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

if [ $1 == call_meth ] ; then
    ##assign meth/unmeth based on given thresholds
    python ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
	   -i $datadir/clin_barocdes_methprobs.csv \
	   -o $datadir/clin_barocdes_methcalls.csv \
	   -m .8 \
	   -u .8
fi



if [ $1 == combine_hiC ] ; then
    ##combine individual hiC bin fastas into a single file
    for i in $clindir/hiC/clusters/*.fasta ;
    do
	name=`basename $i .fasta`
	sed "s/>/>$name./g" $i >> $clindir/hiC/$prefix.hiC.fasta
    done
fi


if [ $1 == mummer_hiC ] ; then
    ##mummer clin sample contigs with bins from phase folks
    mkdir -p $datadir/clin_mummer
    
    nucmer -p $datadir/clin_mummer/asm_hiC $asm $clindir/hiC/$prefix.hiC.fasta
    mummerplot --filter --fat --postscript -p $datadir/clin_mummer/asm_hiC $datadir/clin_mummer/asm_hiC.delta
    mummerplot --filter --fat --png -p $datadir/clin_mummer/asm_hiC $datadir/clin_mummer/asm_hiC.delta
    dnadiff -p $datadir/clin_mummer/asm_hiC $asm $clindir/hiC/$prefix.hiC.fasta
fi
