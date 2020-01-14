#!/bin/bash

##script for everything up to assembly
datadir=/uru/Data/Nanopore/projects/nivar

if [ $1 == r10_untar ] ; then
    prefix=r10_ecoli
    
    tarball=/uru/Data/Nanopore/oxford/191227_ecolik12_r10/191227_ecolik12_r10.tar.gz

    rawdir=$datadir/$prefix/raw
    mkdir -p $rawdir
    
    tar -xzf $tarball -C $rawdir
fi

if [ $1 == r9_untar ] ; then
    prefix=r9_ecoli
    
    tarball=/uru/Data/Nanopore/oxford/191227_ecolik12_r941/191227_ecolik12_r941.tar.gz

    rawdir=$datadir/$prefix/raw
    mkdir -p $rawdir
    
    tar -xzf $tarball -C $rawdir
fi


if [ $1 == r9_call ] ; then
    ##using yfan2012/whack:guppy_bcall
    raw=/media/raw
    guppy_basecaller -i $raw/fast5 -s $raw/called --flowcell FLO-MIN106 --kit SQK-LSK109 -x "cuda:0"
fi


if [ $1 == r10_call ] ; then
    raw=/media/raw
    guppy_basecaller -i $raw/fast5 -s $raw/called --flowcell FLO-MIN110 --kit SQK-LSK109 -x "cuda:0"
fi

##rearranged some of the called/raw data
if [ $1 == gather_fqs ] ; then
    r10=$datadir/r10_ecoli/r10_ecoli.fq
    r9=$datadir/r9_ecoli/r9_ecoli.fq
    
    cat $datadir/r10_ecoli/raw/called/*fastq > $r10
    cat $datadir/r9_ecoli/raw/called/*fastq > $r9

    python ~/Code/utils/fastq_long.py -i $r9 -o $datadir/r9_ecoli/r9_ecoli_3kb.fq -l 3000
    python ~/Code/utils/fastq_long.py -i $r10 -o $datadir/r10_ecoli/r10_ecoli_3kb.fq -l 3000
fi

if [ $1 == assemble ] ; then
    mkdir -p $datadir/assemble

    for i in r9 r10 ;
    do
	mkdir -p $datadir/assemble/${i}_ecoli
	canu \
	    -p ${i}_ecoli -d $datadir/assemble/${i}_ecoli \
	    genomeSize=4.5m \
	    -nanopore-raw $datadir/${i}_ecoli/${i}_ecoli_3kb.fq
    done
fi

if [ $1 == mummer ] ; then
    ref=$datadir/reference/ecoli_k12.fa
    for i in r9 r10 ;
    do
	mkdir -p $datadir/mummer/${i}_ecoli_ref

	raw=$datadir/assemble/${i}_ecoli/${i}_ecoli.contigs.fasta
	
	cp $ref ~/tmp/
	cp $raw ~/tmp/

	nucmer -p ~/tmp/ecoli_${i}_raw_ref ~/tmp/ecoli_k12.fa ~/tmp/${i}_ecoli.contigs.fasta
	mummerplot --filter --fat --png -p ~/tmp/ecoli_${i}_raw_ref ~/tmp/ecoli_${i}_raw_ref.delta
	dnadiff -p ~/tmp/ecoli_${i}_raw_ref ~/tmp/ecoli_k12.fa ~/tmp/${i}_ecoli.contigs.fasta

	cp ~/tmp/ecoli* $datadir/mummer/${i}_ecoli_ref/

	rm ~/tmp/*
	
    done
fi
