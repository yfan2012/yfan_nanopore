#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard/strains
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/allrefs.fa

prefix=220131_mdr_barnyard_
samps='st3294 st3689'

if [ $1 == basecall ] ; then
    mkdir -p $ssddir/called

    for i in $samps ;
    do
	samp=$prefix$i
	mkdir -p $ssddir/called/$samp
	guppy_basecaller \
	    -i $ssddir/raw/$samp/no_sample/*/fast5_pass \
	    -s $ssddir/called/$samp \
	    --recursive \
	    --compress_fastq \
	    --flowcell FLO-FLG001 --kit SQK-RAD004 \
	    --device 'cuda:0'
    done
fi

if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs

    for i in $samps ;
    do
	samp=$prefix$i
	cat $ssddir/called/$samp/pass/*fastq.gz > $datadir/fastqs/$samp.fastq.gz
    done
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/align
    for i in $samps ;
    do
	samp=$prefix$i
	mkdir -p $datadir/align/$samp
	fq=$datadir/fastqs/$samp.fastq.gz

	minimap2 -t 36 -x map-ont $ref $fq \
		 > $datadir/align/$samp.paf
    done
fi

if [ $1 == alignbam ] ; then

    for i in $samps ;
    do
	samp=$prefix$i
	fq=$datadir/fastqs/$samp.fastq.gz

	minimap2 -t 36 -ax map-ont $ref $fq | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/align/$samp.sorted.bam
	samtools index $datadir/align/$samp.sorted.bam
		 
    done
fi


if [ $1 == coverage ] ; then
    for i in $samps ;
    do
	samp=$prefix$i
	bedtools genomecov -d -ibam $datadir/align/$samp.sorted.bam \
		 > $datadir/align/$samp.genomecov
    done
fi


barcodes=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes2.txt
if [ $1 == test ] ; then
    for i in $samps ;
    do
	samp=$prefix$i
	echo $ref
	echo $barcodes
	echo $datadir/megalodon/$samp/per_read_modified_base_calls.txt
	echo $datadir/megalodon/$samp/per_read_modified_base_calls.txt.idx
	echo $datadir/contig_level/$samp.barocdes_methprobs.csv
    done
fi

   
if [ $1 == filter_meth ] ; then
    ##filter out barcode motif related positions in per read meth calls
    mkdir -p $datadir/contig_level
    for i in $samps ;
    do
	samp=$prefix$i
	python3 ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	       -r $ref \
	       -b $barcodes \
	       -m $datadir/megalodon/$samp/per_read_modified_base_calls.txt \
	       -i $datadir/megalodon/$samp/per_read_modified_base_calls.txt.idx \
	       -o $datadir/contig_level/$samp.barocdes_methprobs.csv \
	       -t 12
    done
fi


if [ $1 == call_meth ] ; then
    ##assign meth/unmeth based on given thresholds
    for i in $samps ;
    do
	samp=$prefix$i
	python3 ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
		-i $datadir/contig_level/$samp.barocdes_methprobs.csv \
		-o $datadir/contig_level/$samp.barocdes_methcalls.csv \
		-m .8 \
		-u .8
    done
fi
