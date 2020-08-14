#!/bin/bash

datadir=/uru/Data/Nanopore/projects/mdr
samps="MDRstool_16 MDRstool_19"


if [ $1 == call ] ; then
datadir=~/data/mdr
    for i in $samps ;
    do
	guppy_basecaller \
	    -i $datadir/$i/multiraw \
	    -s $datadir/$i/called \
	    --flowcell FLO-MIN106 --kit SQK-LSK108 \
	    --device 'cuda:0' \
	    --fast5_out
    done
fi

if [ $1 == gather ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/$i/fastqs
	cat $datadir/$i/called/*fastq > $datadir/$i/fastqs/${i}_old.fq
    done
fi

	     
if [ $1 == old_assemble ] ; then
    sizes="150m 100m 500m 10g"
    for i in $samps ;
    do
	for gsize in $sizes ;
	do
	    mkdir -p $datadir/$i/metaflye/$gsize
	    flye \
		--nano-raw $datadir/$i/fastqs/$i.fq \
		-o $datadir/$i/metaflye/$gsize \
		-t 36 \
		-g $gsize \
		--plasmids \
		--meta
	done
    done
fi

if [ $1 == recswab_assemble ] ; then
    gsize=1m
    mkdir -p $datadir/RECswab_1/metaflye/$gsize
    flye \
	--nano-raw $datadir/RECswab_1/fastqs/RECswab_1.fq \
	-o $datadir/RECswab_1/metaflye/$gsize \
	-t 36 \
	-g $gsize \
	--plasmids \
	--meta
fi

	
if [ $1 == rename ] ; then

    for asm in $datadir/*/metaflye/* ;
    do
	echo $asm
	size=`echo $asm | rev | cut -d / -f 1 | rev`
	samp=`echo $asm | rev | cut -d / -f 3 | rev`
	mv $asm/assembly.fasta $asm/${samp}_$size.assembly.fasta
    done
fi


if [ $1 == call_stool16pcr ] ; then
    datadir=~/data/mdr
    guppy_basecaller \
	-i $datadir/multiraw_pcr \
	-s $datadir/called_pcr \
	--flowcell FLO-MIN106 --kit SQK-LSK109 \
	--device 'cuda:0' \
	--fast5_out
fi

if [ $1 == call_stool16native ] ; then
    datadir=~/data/mdr
    guppy_basecaller \
	-i $datadir/200708_mdr_stool16native/fast5 \
	-s $datadir/200708_mdr_stool16native/called \
	--flowcell FLO-MIN106 --kit SQK-LSK109 \
	--device 'cuda:0' \
	--fast5_out
fi

if [ $1 == gather_stool16pcr ] ; then
    cat $datadir/MDRstool_16/called_pcr/*fastq > $datadir/MDRstool_16/fastqs/MDRstool_16_pcr.fq
fi

if [ $1 == gather_stool16native ] ; then
    cat $datadir/MDRstool_16/called_native/*fastq > $datadir/MDRstool_16/fastqs/MDRstool_16_native.fq
fi

##get an estimate of human gut metagenome size
##got zymo human gut ref genomes here: https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip

if [ $1 == cat_refs ] ; then
    cat $datadir/refs/D6331.refseq/genomes/*fasta > $datadir/refs/D6331.refseq/genomes/allrefs.fasta
fi

if [ $1 == check_refsize ] ; then
    python ~/Code/utils/qc/asm_assess.py -i $datadir/refs/D6331.refseq/genomes/allrefs.fasta -p zymo_all
fi

if [ $1 == assemble_stool16 ] ; then
    ##on aws instance
    datadir=~/
    mkdir -p $datadir/metaflye/native_100m
    flye \
	--nano-raw $datadir/fastqs/MDRstool_16_native.fq \
	-o $datadir/metaflye/native_100m \
	-t 36 \
	-g 100m \
	--plasmids \
	--meta

    mkdir -p $datadir/metaflye/pcr_100m
    flye \
	--nano-raw $datadir/fastqs/MDRstool_16_pcr.fq \
	-o $datadir/metaflye/pcr_100m \
	-t 36 \
	-g 100m \
	--plasmids \
	--meta
fi

if [ $1 == rename_stool16 ] ; then
    for i in native_100m pcr_100m ;
    do
	mv $datadir/MDRstool_16/metaflye/$i/assembly.fasta $datadir/MDRstool_16/metaflye/$i/MDRstool_16_$i.assembly.fasta
    done
fi


if [ $1 == assemble ] ; then
    ##try new version of flye that doesn't need gsize param
    for i in native pcr ;
    do
	mkdir -p $datadir/MDRstool_16/metaflye/$i
	
	flye \
	    --nano-raw $datadir/MDRstool_16/fastqs/MDRstool_16_$i.fq \
	    -o $datadir/MDRstool_16/metaflye/$i \
	    -t 36 \
	    --plasmids \
	    --meta

	mv $datadir/MDRstool_16/metaflye/$i/assembly.fasta $datadir/MDRstool_16/metaflye/$i/MDRstool_16_$i.assembly.fasta
    done
fi
