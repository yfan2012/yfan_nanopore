#!/bin/bash

ssddir=~/data/mdr/disco
datadir=/mithril/Data/Nanopore/projects/methbin/disco

#samps='MinION_BF_NAT MinION_CP_NAT MinION_HP_NAT MinION_MH_NAT MinION_NG_NAT MinION_TP_NAT MinION_BA_NAT'
samps=MinION_BA_NAT

ref=$datadir/ref/disco_refs.fasta
if [ $1 == makeref ] ; then
    cat $datadir/ref/*fa > $ref
fi


if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon
    for i in $samps ;
    do
	echo $i
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
    done
fi

megaref=$datadir/ref_meta/meta10.fa
if [ $1 == metagenome ] ; then
    mkdir -p $ssddir/megalodon

    i=MinION_JM3O_NAT
    echo $i
    
    megalodon \
	$ssddir/raw/$i \
	--overwrite \
	--suppress-progress-bars \
	--verbose-read-progress 0 \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	--reference $megaref \
	--outputs per_read_mods \
	--output-directory $ssddir/megalodon/$i \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36     
fi
    
