#!/bin/bash

ssddir=~/data/mdr/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/mdr
dbdir=/atium/Data/ref/mdr_nt
prefix=200708_mdr_stool16native

oldref=$datadir/ref/mdr_refs.fa
ref=$dbdir/mdr_reduced.fa

if [ $1 == megaidx ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt \
	    -o $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt.idx
fi

if [ $1 == megaidx_asm ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $datadir/megalodon/${prefix}_asm/per_read_modified_base_calls.txt \
	    -o $datadir/megalodon/${prefix}_asm/per_read_modified_base_calls.txt.idx
fi

if [ $1 == guppy_call ] ; then
    mkdir -p $ssddir/called

    guppy_basecaller \
	-i $ssddir/raw \
	-s $ssddir/called \
	--recursive \
	--compress_fastq \
	--flowcell FLO-MIN106 --kit SQK-LSK109 \
	--device 'cuda:0'
    
fi

if [ $1 == copy ] ; then
    cp -r $ssddir/called $datadir/
    cp -r $ssddir/megalodon $datadir/
fi

if [ $1 == copy_mega ] ; then
    cp -r $ssddir/megalodon $datadir/
fi

if [ $1 == barcode ] ; then
    mkdir -p $datadir/barcode
    mkdir -p $datadir/barcode/$prefix
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
           -o $datadir/barcode/$prefix/${prefix}_barcodes.txt \
           -t 12 ;} &> $datadir/barcode/$prefix/${prefix}_time.txt
fi

if [ $1 == barcode_asm ] ; then
    mkdir -p $datadir/barcode
    mkdir -p $datadir/barcode/${prefix}_asm
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/${prefix}_asm/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/${prefix}_asm/per_read_modified_base_calls.txt.idx \
           -r $datadir/medaka/consensus.fasta \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
           -o $datadir/barcode/${prefix}_asm/${prefix}_barcodes50.txt \
           -t 12 ;} &> $datadir/barcode/${prefix}_asm/${prefix}_time.txt
fi


if [ $1 == align ] ; then
    mkdir -p $datadir/align
    
    fq=$datadir/fastqs/$prefix.fq.gz
    minimap2 -t 36 -x map-ont $ref $fq \
	     > $datadir/align/$prefix.paf
fi

if [ $1 == align_asm ] ; then
    mkdir -p $datadir/align
    
    fq=$datadir/fastqs/$prefix.fq.gz
    minimap2 -t 36 -x map-ont $datadir/medaka/consensus.fasta $fq \
	     > $datadir/align/${prefix}_asm.paf
fi


if [ $1 == motifcounts ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -r $ref \
	   -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
	   -a $datadir/align/$prefix.paf \
	   -m $datadir/barcode/$prefix/${prefix}_barcodes.txt \
	   -o $datadir/barcode/$prefix/${prefix}_barcodes_motifcounts.txt \
	   -q 40 \
	   -v
fi

if [ $1 == motifcounts_asm ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -r $datadir/medaka/consensus.fasta \
	   -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
	   -a $datadir/align/${prefix}_asm.paf \
	   -m $datadir/barcode/${prefix}_asm/${prefix}_barcodes50.txt \
	   -o $datadir/barcode/${prefix}_asm/${prefix}_barcodes50_motifcounts.txt \
	   -q 40 \
	   -v
fi


if [ $1 == alignbam ] ; then
    fq=$datadir/fastqs/$prefix.fq.gz
    minimap2 -ax map-ont $ref $fq | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/align/$prefix.sorted.bam
    samtools index $datadir/align/$prefix.sorted.bam
fi


if [ $1 == barcode50 ] ; then
    mkdir -p $datadir/barcode
    mkdir -p $datadir/barcode/$prefix
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
           -o $datadir/barcode/$prefix/${prefix}_barcodes_50.txt \
           -t 12 ;} &> $datadir/barcode/$prefix/${prefix}_time_50.txt
fi

if [ $1 == motifcounts50 ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -r $ref \
	   -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
	   -a $datadir/align/$prefix.paf \
	   -m $datadir/barcode/$prefix/${prefix}_barcodes_50.txt \
	   -o $datadir/barcode/$prefix/${prefix}_barcodes_motifcounts_50.txt \
	   -q 40 \
	   -v
fi
