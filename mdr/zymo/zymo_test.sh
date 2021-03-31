#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/raw
datadir=~/data/mdr/zymo
prefix=20190809_zymo_control

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw/$prefix
    tar -xzf $rawdir/$prefix.tar.gz -C $datadir/raw/$prefix
fi

##ref from https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
ref=/uru/Data/Nanopore/projects/read_class/ref/zymo_all.fa
if [ $1 == makeref ] ; then
    cat /uru/Data/Nanopore/projects/read_class/ref/D6322.refseq/Genomes/*fasta > $ref
fi
    
if [ $1 == megalodon ] ; then
    mkdir -p $datadir/megalodon
    mkdir -p $datadir/megalodon/$prefix
    
    megalodon \
        $datadir/raw/$prefix/$prefix/fast5 \
        --overwrite \
        --guppy-server-path "/usr/bin/guppy_basecall_server" \
        --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
        --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
        --reference $ref \
        --outputs per_read_mods \
        --output-directory $datadir/megalodon/$prefix \
        --write-mods-text \
        --devices "cuda:0" \
        --processes 36
    rm $datadir/megalodon/$prefix/per_read_modified_base_calls.db
fi

if [ $1 == index ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
	    -o $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx
fi

	    
if [ $1 == call ] ; then
    mkdir -p $datadir/barcode
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/qc/barcodes.txt \
           -o $datadir/barcode/zymo_barcodes.txt \
	   -n 100000000000 \
           -t 36 ;} &> $datadir/barcode/zymo_time.txt
fi
