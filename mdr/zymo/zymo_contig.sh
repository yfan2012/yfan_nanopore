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
	--meta
    mv $datadir/flye/$prefix/assembly.fasta $datadir/flye/$prefix/$prefix.fasta
    ##was rerun with v2.9
fi

fq=$datadir/fastq/$prefix/$prefix.fq.gz

if [ $1 == racon_align ] ; then
    mkdir -p $datadir/racon

    minimap2 -t 36 -ax map-ont $datadir/flye/$prefix/$prefix.fasta $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $datadir/racon/$prefix.sorted.bam
    samtools index $datadir/racon/$prefix.sorted.bam
fi

if [ $1 == racon ] ; then		   
    samtools view -@ 36 $datadir/racon/$prefix.sorted.bam > $datadir/racon/$prefix.sorted.sam
    racon -m 8 -x -6 -g -8 -w 500 -t 54\
	  $fq \
	  $datadir/racon/$prefix.sorted.sam \
	  $datadir/flye/$prefix/$prefix.fasta > $datadir/racon/$prefix.racon.contigs.fasta
fi


if [ $1 == medaka ] ; then
    mkdir -p $datadir/medaka

    medaka_consensus \
	-i $fq \
	-d $datadir/racon/$prefix.racon.contigs.fasta \
	-o $datadir/medaka \
	-t 54 \
	-m r941_min_high_g360
fi


if [ $1 == untar ] ; then
    mkdir -p $ssddir/raw/$prefix
    tar -xzf $rawdir/$prefix.tar.gz -C $ssddir/raw/$prefix
fi


if [ $1 == megalodon_polished ] ; then
    mkdir -p $ssddir/megalodon

    megalodon \
        $ssddir/raw/$prefix/$prefix/fast5 \
        --overwrite \
        --guppy-server-path "/usr/bin/guppy_basecall_server" \
        --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
        --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
        --reference $datadir/medaka/consensus.fasta \
        --outputs per_read_mods \
        --output-directory $ssddir/megalodon/${prefix}_polished \
        --write-mods-text \
        --devices "cuda:0" \
        --processes 36
fi

if [ $1 == move_polished ] ; then
    mkdir -p $datadir/megalodon/${prefix}_polished

    mv $ssddir/megalodon/${prefix}_polished/* $datadir/megalodon/${prefix}_polished/
fi

if [ $1 == megaidx_polished ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt \
	    -o $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt.idx
fi

if [ $1 == barcode_polished ] ; then
    mkdir -p $datadir/barcode/${prefix}_polished
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt.idx \
           -r $datadir/medaka/consensus.fasta \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
           -o $datadir/barcode/${prefix}_polished/${prefix}_contigs_barcodes15.txt \
           -t 12 ;} &> $datadir/barcode/${prefix}_polished/${prefix}_contigs_barcodes15_time.txt
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

if [ $1 == move_meta ] ; then
    mkdir -p $datadir/megalodon/${prefix}_meta

    mv $ssddir/megalodon/${prefix}_meta/* $datadir/megalodon/${prefix}_meta/
fi


if [ $1 == megaidx ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt \
	    -o $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt.idx
fi

if [ $1 == barcode_common ] ; then
    mkdir -p $datadir/barcode/${prefix}_meta
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt.idx \
           -r $asm \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
           -o $datadir/barcode/${prefix}_meta/${prefix}_contigs_barcodes15.txt \
           -t 36 ;} &> $datadir/barcode/${prefix}_meta/${prefix}_contigs_barcodes15_time.txt
fi
    
