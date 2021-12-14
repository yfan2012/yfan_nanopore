#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/zymo/raw
ssddir=~/data/mdr/zymo
prefix=20190809_zymo_control
projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/etc

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


if [ $1 == move_stuff ] ; then
    cp -r $ssddir/megalodon/CG_model $datadir/
fi
	    
if [ $1 == extract_motifs ] ; then
    python ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $ref \
	   -b ~/Code/yfan_nanopore/mdr/zymo/barcodes_zymo_curated.txt \
	   -m $datadir/CG_model/per_read_modified_base_calls.txt \
	   -i $datadir/CG_model/per_read_modified_base_calls.txt.idx \
	   -o $datadir/CG_model/curate_filter.csv \
	   -t 12
fi
    
if [ $1 == filter_motifs ] ; then
    python ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
	   -i $datadir/CG_model/curate_filter.csv \
	   -o $datadir/CG_model/curate_calls.csv \
	   -m .8 \
	   -u .8
fi


if [ $1 == check_motifs ] ; then

    while read p ;
    do
	label=`echo $p | cut -d ' ' -f 1` 
	motifs=`echo $p | cut -d ' ' -f 2`
	positions=`echo $p | cut -d ' ' -f 3`
	echo $label
	python ~/Code/yfan_meth/utils/bismark_motif_confirm.py \
	       -c $projdir/zymo/truth/bisulfite/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
	       -r $ref \
	       -m $motifs \
	       -p $positions
    done < bisulfinfo.txt
fi


fq=$projdir/zymo/fastq/$prefix/$prefix.fq.gz
if [ $1 == nanopolishidx ] ; then
    nanopolish index \
	       -d $ssddir/raw/$prefix/$prefix/fast5 \
	       -s $projdir/zymo/called/$prefix/sequencing_summary.txt \
	       $fq
fi

if [ $1 == align ] ; then
    minimap2 -t 36 -ax map-ont $ref $fq | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/CG_model/$prefix.sorted.bam
    samtools index $datadir/CG_model/$prefix.sorted.bam
fi
	
if [ $1 == nanopolish ] ; then
    nanopolish call-methylation \
	       -r $fq \
	       -b $datadir/CG_model/$prefix.sorted.bam \
	       -g $ref \
	       -q cpg \
	       -t 36 > $datadir/CG_model/nanopolish_methcalls.tsv
fi
    
if [ $1 == nanopolishfreq ] ; then
    ~/software/nanopolish/scripts/calculate_methylation_frequency.py -s \
	$datadir/CG_model/nanopolish_methcalls.tsv \
	> $datadir/CG_model/nanopolish_methcalls_freq.tsv
fi


if [ $1 == nanopolishfreq_nosplit ] ; then
    ~/software/nanopolish/scripts/calculate_methylation_frequency.py \
	$datadir/CG_model/nanopolish_methcalls.tsv \
	> $datadir/CG_model/nanopolish_methcalls_freq_nosplit.tsv
fi
