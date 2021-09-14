#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/allrefs.fa

##prefix=210908_mdr_barnyard_mix1
##prefix=210908_mdr_barnyard_mix2
##prefix=210912_mdr_barnyard_mix3
prefix=210912_mdr_barnyard_mix4

if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon
    mkdir -p $ssddir/megalodon/$prefix

    megalodon \
	$ssddir/raw/$prefix/no_sample/*/fast5_pass \
	--overwrite \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	--reference $ref \
	--outputs per_read_mods \
	--output-directory $ssddir/megalodon/$prefix \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36
    rm $ssddir/megalodon/$prefix/per_read_modified_base_calls.db
fi

if [ $1 == move ] ; then
    mkdir -p $datadir/megalodon
    mv $ssddir/megalodon/$prefix $datadir/megalodon/
fi
