#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/speciescheck_refs.fa

prefix=210922_mdr_barnyard_speciescheck

if [ $1 == addstaph ] ; then
    ##grab staph read names from the no plasmid mix
    awk '$6=="NC_007795.1" {print $1}' $datadir/align/210908_mdr_barnyard_mix1.paf > $ssddir/raw/$prefix/staph_readnames.txt

    fast5_subset \
	-i $ssddir/raw/210908_mdr_barnyard_mix1/no_sample/*/fast5_pass \
	-s $ssddir/raw/$prefix/no_sample/*/fast5_pass/staph_reads \
	-l $ssddir/raw/$prefix/staph_readnames.txt \
	-f staph_mix1 
	
fi

if [ $1 == basecall ] ; then
    mkdir -p $ssddir/called
    mkdir -p $ssddir/called/$prefix
    
    guppy_basecaller \
	-i $ssddir/raw/$prefix \
	-s $ssddir/called/$prefix \
	--recursive \
	--compress_fastq \
	--flowcell FLO-FLG001 --kit SQK-RAD004 \
	--device 'cuda:0'
fi

if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs

    cat $ssddir/called/$prefix/pass/*fastq.gz > $datadir/fastqs/$prefix.fastq.gz
fi


if [ $1 == align ] ; then
    mkdir -p $datadir/align

    fq=$datadir/fastqs/$prefix.fastq.gz
    
    minimap2 -t 36 -x map-ont $ref $fq \
	     > $datadir/align/$prefix.paf
	
fi

if [ $1 == alignbam ] ; then

    fq=$datadir/fastqs/$prefix.fastq.gz
    
    minimap2 -t 36 -ax map-ont $ref $fq | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/align/$prefix.sorted.bam
    samtools index $datadir/align/$prefix.sorted.bam
    
fi

	
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

if [ $1 == megaidx ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	     -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
	     -o $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx
fi

if [ $1 == barcode ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode.py \
            -m $datadir/megalodon/$prefix/per_read_modified_base_calls.txt \
            -i $datadir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
            -r $ref \
            -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
            -o $datadir/barcode/${prefix}_barcodes.txt \
            -t 12
fi


if [ $1 == motifcounts ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	    -r $ref \
	    -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
	    -a $datadir/align/$prefix.paf \
	    -m $datadir/barcode/${prefix}_barcodes.txt \
	    -o $datadir/barcode/${prefix}_barcodes_motifcounts.txt \
	    -q 40
fi
