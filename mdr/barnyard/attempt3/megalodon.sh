#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/allrefs.fa

prefix=211005_mdr_barnyard_mix

if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon

    for i in 5 6 ;
    do
	samp=${prefix}$i
	mkdir -p $ssddir/megalodon/$samp
	
	megalodon \
	    $ssddir/raw/$samp/no_sample/*/fast5_pass \
	    --overwrite \
	    --guppy-server-path "/usr/bin/guppy_basecall_server" \
	    --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	    --reference $ref \
	    --outputs per_read_mods \
	    --output-directory $ssddir/megalodon/$samp \
	    --write-mods-text \
	    --devices "cuda:0" \
	    --processes 36
	rm $ssddir/megalodon/$samp/per_read_modified_base_calls.db
    done
fi

if [ $1 == move ] ; then
    mkdir -p $datadir/megalodon
    for i in 5 6 ;
    do
	samp=${prefix}$i
	mv $ssddir/megalodon/$samp $datadir/megalodon/
    done
fi
