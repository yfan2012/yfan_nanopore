#!/bin/bash

projdir=/mithril/Data/Nanopore/projects/methbin/zymo
datadir=$projdir/contig_agg
prefix=20190809_zymo_control
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == test ] ; then
    mkdir -p $datadir
    mkdir -p $datadir/test

    time ( python ~/Code/yfan_meth/utils/megalodon_agg.py \
	   -m $projdir/megalodon/$prefix/per_read_modified_base_calls.txt \
	   -i $projdir/megalodon/$prefix/per_read_modified_base_calls.txt.small.idx \
	   -r $ref \
	   -o $datadir/test/test.txt \
	   -t 12 \
	   -v )
fi

if [ $1 == agg ] ; then
    mkdir -p $datadir
    mkdir -p $datadir/$prefix

    python ~/Code/yfan_meth/utils/megalodon_agg.py \
	   -m $projdir/megalodon/$prefix/per_read_modified_base_calls.txt \
	   -i $projdir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
	   -r $ref \
	   -o $datadir/$prefix/$prefix.meth_report.txt \
	   -t 12 \
	   -v
fi

if [ $1 == mega_agg ] ; then
    ssddir=~/data/mdr/zymo
    megalodon \
	$ssddir/raw/$prefix/$prefix/fast5 \
	--overwrite \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	--reference $ref \
	--outputs mods mod_mappings \
	--output-directory $ssddir/megalodon/$prefix \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36
fi



if [ $1 == findmotifs ] ; then
    mkdir -p $datadir/$prefix

    while read p; do
        label=`echo $p | cut -d ' ' -f 2`
        chrom=`echo $p | cut -d ' ' -f 1`
	python ~/Code/yfan_meth/utils/bismark_motif_finder.py \
	       -c $datadir/$prefix/$prefix.meth_report.txt \
	       -r $ref \
	       -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
	       -p .5 \
	       -l 8 \
	       -m 10 \
	       -s $chrom >> $datadir/$prefix/${prefix}_barcodes50.csv
	done < ./truth/chrlist_withlabels.txt
fi


if [ $1 == extract_motifs ] ; then
    ##extract meth info at positions relelvant to the motifs in the barcode
    
    time (python3 ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $ref \
	   -b ~/Code/yfan_nanopore/mdr/zymo/barcodes_zymo_curated.txt \
	   -m $projdir/megalodon/$prefix/per_read_modified_base_calls.txt \
	   -i $projdir/megalodon/$prefix/per_read_modified_base_calls.txt.idx \
	   -o $datadir/$prefix/$prefix.curate_extract.csv \
	   -t 12 )
fi

    
if [ $1 == filter_motifs ] ; then
    time( python ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
		 -i $datadir/$prefix/$prefix.curate_extract.csv \
		 -o $datadir/$prefix/$prefix.curate_filter.csv \
		 -m .8 \
		 -u .8)
fi

fq=$projdir/fastq/$prefix/$prefix.fq.gz
if [ $1 == align ] ; then
    mkdir -p $datadir/align

    minimap2 -t 36 -ax map-ont $ref $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $datadir/align/${prefix}.sorted.bam
    samtools index $datadir/align/${prefix}.sorted.bam
fi

if [ $1 == coverage ] ; then
    bedtools genomecov -d \
	     -ibam $datadir/align/${prefix}.sorted.bam \
	     > $datadir/align/${prefix}.sorted.cov
fi
