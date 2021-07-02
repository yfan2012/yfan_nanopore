#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/zymo/raw
ssddir=~/data/mdr/zymo
prefix=20190809_zymo_control
datadir=/mithril/Data/Nanopore/projects/methbin


if [ $1 == untar ] ; then
    mkdir -p $datadir/raw/$prefix
    tar -xzf $rawdir/$prefix.tar.gz -C $ssddir/raw/$prefix
fi

##ref from https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa
if [ $1 == makeref ] ; then
    cat /uru/Data/Nanopore/projects/read_class/ref/D6322.refseq/Genomes/*fasta > $ref
fi
    
if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon
    mkdir -p $ssddir/megalodon/$prefix
    
    megalodon \
        $ssddir/raw/$prefix/$prefix/fast5 \
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


if [ $1 == index ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt \
	    -o $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt.idx
fi

	    
if [ $1 == call ] ; then
    mkdir -p $ssddir/barcode
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $ssddir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/zymo/zymo_barcodes.txt \
           -o $ssddir/barcode/${prefix}_barcodes.txt \
	   -n 100000000000 \
           -t 54 ;} &> $ssddir/barcode/${prefix}_time.txt
fi


if [ $1 == basecall ] ; then
    mkdir -p $ssddir/called

    guppy_basecaller \
	-i $ssddir/raw/$prefix/$prefix/fast5 \
	-s $ssddir/called/$prefix \
	--compress_fastq \
	--flowcell FLO-MIN106 --kit SQK-LSK109 \
	--device 'cuda:0'
fi


if [ $1 == gather ] ; then
    mkdir -p $ssddir/fastq/$prefix

    cat $ssddir/called/$prefix/pass/*fastq.gz > $ssddir/fastq/$prefix/$prefix.fq.gz
fi


if [ $1 == align ] ; then
    mkdir -p $ssddir/align
    
    fq=$ssddir/fastq/$prefix/$prefix.fq.gz
    minimap2 -t 36 -x map-ont $ref $fq \
	     > $ssddir/align/$prefix.paf
fi


if [ $1 == filter ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -a $ssddir/align/$prefix.paf \
	   -r $ref \
	   -o $ssddir/barcode/${prefix}_barcodes_filtered.txt \
	   -m $ssddir/barcode/${prefix}_barcodes.txt \
	   -b ~/Code/yfan_nanopore/mdr/zymo/zymo_barcodes.txt \
	   -q 40 \
	   -l 5000 \
	   -n 5 \
	   -v
fi

if [ $1 == call_moved ] ; then
    mkdir -p $datadir/barcode
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/zymo/zymo_barcodes.txt \
           -o $datadir/zymo/barcode/${prefix}_barcodes.txt \
           -t 36 ;} &> $datadir/zymo/barcode/${prefix}_time.txt
fi
    
if [ $1 == filtercount ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -a $datadir/zymo/align/$prefix.paf \
	   -r $ref \
	   -o $datadir/zymo/barcode/${prefix}_motifcounts.txt \
	   -m $datadir/zymo/barcode/${prefix}_barcodes.txt \
	   -b ~/Code/yfan_nanopore/mdr/zymo/zymo_barcodes.txt \
	   -q 40 \
	   -v
fi

if [ $1 == callcommon ] ; then
    mkdir -p $datadir/zymo/barcode
    #for i in 10 15 20 ;
    for i in 20 ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt \
               -i $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/rebase/barcodes${i}.txt \
               -o $datadir/zymo/barcode/${prefix}_barcodes${i}.txt \
	       -n 100000000000 \
               -t 12 ;} &> $datadir/zymo/barcode/${prefix}_time${i}.txt
    done
fi

if [ $1 == filtcommon ] ; then
    #for i in 10 15 20 ;
    for i in 20 ;
    do
	python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	       -a $datadir/zymo/align/$prefix.paf \
	       -r $ref \
	       -o $datadir/zymo/barcode/${prefix}_motifcounts${i}.txt \
	       -m $datadir/zymo/barcode/${prefix}_barcodes${i}.txt \
	       -b ~/Code/yfan_nanopore/mdr/rebase/barcodes${i}.txt \
	       -q 40 \
	       -v
    done
fi


if [ $1 == contig_common_test ] ; then
    i=15
    head -100 $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
	 > $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt.small.idx
    
    { time python ~/Code/yfan_meth/utils/megalodon_contig.py \
           -m $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt \
           -i $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt.small.idx \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes${i}.txt \
           -o $datadir/zymo/contig/${prefix}_contig${i}.txt \
	   -v $datadir/zymo/contig/${prefix}_contig${i}_cov.txt \
	   -c 1 \
	   -a 0 \
           -t 12 ;} &> $datadir/zymo/contig/${prefix}_time${i}.small.txt
fi
    

if [ $1 == contig_common ] ; then
    mkdir -p $datadir/zymo/contig
    for i in 15 ;
    do
	{ time python ~/Code/yfan_meth/utils/megalodon_contig.py \
               -m $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt \
               -i $datadir/zymo/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/rebase/barcodes${i}.txt \
               -o $datadir/zymo/contig/${prefix}_contig${i}.txt \
	       -v $datadir/zymo/contig/${prefix}_contig${i}_cov.txt \
	       -c 1 \
	       -a 0 \
               -t 12 ;} &> $datadir/zymo/contig/${prefix}_time${i}.txt
    done
fi
