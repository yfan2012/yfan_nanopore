#!/bin/bash

prefix=200708_mdr_stool16native
projdir=/mithril/Data/Nanopore/projects/methbin
clindir=$projdir/mdr
datadir=$projdir/paperfigs/contig_level

barcodes=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes.txt
barcodes2=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes2.txt
barcodes3=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes3.txt
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
    python3 ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
	   -i $datadir/clin_barocdes_methprobs.csv \
	   -o $datadir/clin_barocdes_methcalls.csv \
	   -m .8 \
	   -u .8
fi

if [ $1 == meth_perf_idx ] ; then
    ##aggregate meth calls from the alignment filtered reads
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	   -i $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt \
	   -o $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt.idx

    ##filter out barcode motif related positions in per read meth calls
    python3 ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $asm \
	   -b $barcodes \
	   -m $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt \
	   -i $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt.idx \
	   -o $datadir/clin_barocdes_methprobs.perf.csv \
	   -t 12
fi

if [ $1 == meth_perf_filter ] ; then
    ##assign meth/unmeth based on given thresholds
    python3 ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
	   -i $datadir/clin_barocdes_methprobs.perf.csv \
	   -o $datadir/clin_barocdes_methcalls.perf.csv \
	   -m .8 \
	   -u .8
fi


if [ $1 == meth_perf2 ] ; then
    ##filter out barcode motif related positions in per read meth calls
    python3 ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $asm \
	   -b $barcodes2 \
	   -m $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt \
	   -i $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt.idx \
	   -o $datadir/clin_barocdes_methprobs.perf2.csv \
	   -t 12

fi

if [ $1 == meth_perf_filter2 ] ; then
    ##assign meth/unmeth based on given thresholds
    python3 ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
	   -i $datadir/clin_barocdes_methprobs.perf2.csv \
	   -o $datadir/clin_barocdes_methcalls.perf2.csv \
	   -m .8 \
	   -u .8
fi

if [ $1 == meth_perf3 ] ; then
    ##filter out barcode motif related positions in per read meth calls
    python3 ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $asm \
	   -b $barcodes3 \
	   -m $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt \
	   -i $clindir/megalodon/${prefix}_perf/per_read_modified_base_calls.txt.idx \
	   -o $datadir/clin_barocdes_methprobs.perf3.csv \
	   -t 12

fi

if [ $1 == meth_perf_filter3 ] ; then
    ##assign meth/unmeth based on given thresholds
    python3 ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
	   -i $datadir/clin_barocdes_methprobs.perf3.csv \
	   -o $datadir/clin_barocdes_methcalls.perf3.csv \
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

