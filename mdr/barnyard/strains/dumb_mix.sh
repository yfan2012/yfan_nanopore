#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard/strains
ssddir=~/data/mdr/barnyard
ref=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/allrefs.fa

mixes='both_plas ecoli_plas staph_plas'

if [ $1 == sepraw ] ; then
    awk '$6!="PRW62" {print $1}' $datadir/align/220131_mdr_barnyard_st3294.paf \
	> $datadir/staph_chr_list.txt
    awk '$6=="PRW62" {print $1}' $datadir/align/220131_mdr_barnyard_st3294.paf \
	> $datadir/staph_plas_list.txt
    awk '$6!="PRW62" {print $1}' $datadir/align/220131_mdr_barnyard_st3689.paf \
	> $datadir/ecoli_chr_list.txt
    awk '$6=="PRW62" {print $1}' $datadir/align/220131_mdr_barnyard_st3689.paf \
	> $datadir/ecoli_plas_list.txt
fi

if [ $1 == copy ] ; then
    ##from coverage analysis:
    ##2.02x more staph cov
    ##13.42 more ecoli plas
    
    mkdir -p $ssddir/raw/both_plas
    mkdir -p $ssddir/raw/ecoli_plas
    mkdir -p $ssddir/raw/staph_plas


    ##copy all ecoli chr reads to all samps
    for samp in $mixes ;
    do
	fast5_subset \
	    -i $ssddir/raw/220131_mdr_barnyard_st3689/no_sample/*/fast5_pass \
	    -s $ssddir/raw/$samp \
	    -f echr \
	    -l $datadir/ecoli_chr_list.txt
    done

    ##copy all staph plas reads to the right places
    for samp in staph_plas both_plas ;
    do
	fast5_subset \
	    -i $ssddir/raw/220131_mdr_barnyard_st3294/no_sample/*/fast5_pass \
	    -s $ssddir/raw/$samp \
	    -f splas \
	    -l $datadir/staph_plas_list.txt
    done

    ##copy just the first half of staph chr reads to all samps
    numreads=`wc -l $datadir/staph_chr_list.txt | awk '{print $1}'`
    takenum=`awk -v numreads=$numreads 'BEGIN { print numreads/2.02 }'`
    head -n $takenum $datadir/staph_chr_list.txt > $datadir/staph_chr_sublist.txt
    for samp in $mixes ;
    do
	fast5_subset \
	    -i $ssddir/raw/220131_mdr_barnyard_st3294/no_sample/*/fast5_pass \
	    -s $ssddir/raw/$samp \
	    -f schr \
	    -l $datadir/staph_chr_sublist.txt
    done

    numreads=`wc -l $datadir/ecoli_plas_list.txt | awk '{print $1}'`
    takenum=`awk -v numreads=$numreads 'BEGIN { print int(numreads/13.42) }'`
    head -n $takenum $datadir/ecoli_plas_list.txt > $datadir/ecoli_plas_sublist.txt
    for samp in ecoli_plas both_plas ;
    do
	fast5_subset \
	    -i $ssddir/raw/220131_mdr_barnyard_st3689/no_sample/*/fast5_pass \
	    -s $ssddir/raw/$samp \
	    -f eplas \
	    -l $datadir/ecoli_plas_sublist.txt
    done
fi
    

if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon

    for samp in $mixes ;
    do
	mkdir -p $ssddir/megalodon/$samp
	
	megalodon \
	    $ssddir/raw/$samp \
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

if [ $1 == movemega ] ; then
    mkdir -p $datadir/megalodon
    for samp in $mixes ;
    do
	mv $ssddir/megalodon/$samp $datadir/megalodon/
    done
fi


if [ $1 == megaidx ] ; then
    for samp in $mixes ;
    do
	python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
		-i $datadir/megalodon/$samp/per_read_modified_base_calls.txt \
		-o $datadir/megalodon/$samp/per_read_modified_base_calls.txt.idx
    done
fi


if [ $1 == basecall ] ; then
    for samp in $mixes ;
    do
	mkdir -p $ssddir/called/$samp
	guppy_basecaller \
	    -i $ssddir/raw/$samp \
	    -s $ssddir/called/$samp \
	    --recursive \
	    --compress_fastq \
	    --flowcell FLO-FLG001 --kit SQK-RAD004 \
	    --device 'cuda:0'

    done
fi

if [ $1 == mvbcall ] ; then
    for samp in $mixes ;
    do
	mv $ssddir/called/$samp $datadir/called/
    done
fi


if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs

    for samp in $mixes ;
    do
	cat $datadir/called/$samp/pass/*fastq.gz > $datadir/fastqs/$samp.fastq.gz
    done
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/align
    for samp in $mixes ;
    do
	mkdir -p $datadir/align/$samp
	fq=$datadir/fastqs/$samp.fastq.gz

	minimap2 -t 36 -x map-ont $ref $fq \
		 > $datadir/align/$samp.paf
    done
fi

if [ $1 == alignbam ] ; then

    for samp in $mixes ;
    do
	fq=$datadir/fastqs/$samp.fastq.gz

	minimap2 -t 36 -ax map-ont $ref $fq | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/align/$samp.sorted.bam
	samtools index $datadir/align/$samp.sorted.bam
		 
    done
fi


if [ $1 == genomecov ] ; then
    for samp in $mixes ;
    do
	bedtools genomecov -d -ibam $datadir/align/$samp.sorted.bam \
		 > $datadir/align/$samp.genomecov
    done
fi


barcodes=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes2.txt
if [ $1 == test ] ; then
    for samp in $mixes ;
    do
	echo $ref
	echo $barcodes
	echo $datadir/megalodon/$samp/per_read_modified_base_calls.txt
	echo $datadir/megalodon/$samp/per_read_modified_base_calls.txt.idx
	echo $datadir/contig_level/$samp.barocdes_methprobs.csv
    done
fi

   
if [ $1 == filter_meth ] ; then
    ##filter out barcode motif related positions in per read meth calls
    mkdir -p $datadir/contig_level
    for samp in $mixes ;
    do
	python3 ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	       -r $ref \
	       -b $barcodes \
	       -m $datadir/megalodon/$samp/per_read_modified_base_calls.txt \
	       -i $datadir/megalodon/$samp/per_read_modified_base_calls.txt.idx \
	       -o $datadir/contig_level/$samp.barocdes_methprobs.csv \
	       -t 12
    done
fi


if [ $1 == call_meth ] ; then
    ##assign meth/unmeth based on given thresholds
    for samp in $mixes ;
    do
	python3 ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
		-i $datadir/contig_level/$samp.barocdes_methprobs.csv \
		-o $datadir/contig_level/$samp.barocdes_methcalls.csv \
		-m .8 \
		-u .8
    done
fi

