#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/zymo/raw
ssddir=~/data/mdr/zymo
prefix=20190809_zymo_control
datadir=/mithril/Data/Nanopore/projects/methbin/zymo


##ref from https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == assemble ] ; then
    mkdir -p $datadir/flye

    flye \
	--nano-raw $datadir/fastq/$prefix/$prefix.fq.gz \
	-o $datadir/flye/$prefix \
	-t 36 \
	-g 100m \
	--plasmids \
	--meta
    mv $datadir/flye/$prefix/assembly.fasta $datadir/flye/$prefix/$prefix.fasta
fi


if [ $1 == untar ] ; then
    mkdir -p $ssddir/raw/$prefix
    tar -xzf $rawdir/$prefix.tar.gz -C $ssddir/raw/$prefix
fi


asm=$datadir/flye/$prefix/$prefix.fasta
if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon

    megalodon \
        $ssddir/raw/$prefix/$prefix/fast5 \
        --overwrite \
        --guppy-server-path "/usr/bin/guppy_basecall_server" \
        --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
        --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
        --reference $asm \
        --outputs per_read_mods \
        --output-directory $ssddir/megalodon/${prefix}_meta \
        --write-mods-text \
        --devices "cuda:0" \
        --processes 36
fi
