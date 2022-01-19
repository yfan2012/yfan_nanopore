#!/bin/bash

ssddir=~/data/mdr/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/mdr
i=200708_mdr_stool16native


ref=$datadir/ref/mdr_refs.fa
if [ $1 == makeref ] ; then
    ##probs not strictly needed idk
    gunzip $datadir/ref/mdr_refs.fa.gz
fi


if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon

    megalodon \
	$ssddir/raw/$i \
	--overwrite \
	--suppress-progress-bars \
	--verbose-read-progress 0 \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	--reference $ref \
	--outputs per_read_mods \
	--output-directory $ssddir/megalodon/$i \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36     
    
fi


asm=$datadir/medaka/consensus.fasta
if [ $1 == megalodon_asm ] ; then
    ##run megalodon with the metagenomic asm
    mkdir -p $ssddir/megalodon

    megalodon \
	$ssddir/raw/$i \
	--overwrite \
	--suppress-progress-bars \
	--verbose-read-progress 0 \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	--reference $asm \
	--outputs per_read_mods \
	--output-directory $ssddir/megalodon/${i}_asm \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36     
fi

if [ $1 == move_asm ] ; then
    rm $ssddir/megalodon/${i}_asm/per_read_modified_base_calls.db

    mkdir -p $datadir/megalodon/${i}_asm
    mv $ssddir/megalodon/${i}_asm/* $datadir/megalodon/${i}_asm/
fi




if [ $1 == find_perf_fast5 ] ; then
    fast5_subset \
	-i $ssddir/raw/$i \
	-s $ssddir/raw/${i}_perf \
	-l $datadir/align/${i}_perf_readids.txt \
	--recursive
fi


if [ $1 == megalodon_perf ] ; then
    ##run megalodon with the asm and reads that only align uniquely with perf mapq
    megalodon \
	$ssddir/raw/${i}_perf \
	--overwrite \
	--suppress-progress-bars \
	--verbose-read-progress 0 \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	--reference $asm \
	--outputs per_read_mods \
	--output-directory $ssddir/megalodon/${i}_perf \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36
    rm $ssddir/megalodon/${i}_perf/per_read_modified_base_calls.db
    

    mkdir -p $datadir/megalodon/${i}_perf
    mv $ssddir/megalodon/${i}_perf/* $datadir/megalodon/${i}_perf/
fi
    
