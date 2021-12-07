#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/zymo/raw
ssddir=~/data/mdr/zymo
prefix=20190809_zymo_control
datadir=/mithril/Data/Nanopore/projects/methbin

##ref from https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa
    
if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon
    mkdir -p $ssddir/megalodon/$prefix
    
    megalodon \
        $ssddir/raw/$prefix/$prefix/fast5 \
        --overwrite \
        --guppy-server-path "/usr/bin/guppy_basecall_server" \
        --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
        --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg \
        --reference $ref \
	--mod-motif m CG 0 \
        --outputs per_read_mods mods mod_mappings \
        --output-directory $ssddir/megalodon/CG_model \
        --write-mods-text \
        --devices "cuda:0" \
        --processes 36
    rm $ssddir/megalodon/CG_model/per_read_modified_base_calls.db
fi


if [ $1 == index ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $ssddir/megalodon/CG_model/per_read_modified_base_calls.txt \
	    -o $ssddir/megalodon/CG_model/per_read_modified_base_calls.txt.idx
fi

	    
