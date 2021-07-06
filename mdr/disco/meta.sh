#!/bin/bash

ssddir=~/data/mdr/disco
datadir=/mithril/Data/Nanopore/projects/methbin/disco
samps=MinION_JM3O_NAT

ref=$datadir/ref_meta/meta10.fa

if [ $1 == makeref ] ; then
    cat $datadir/ref/*fa > $ref
fi

if [ $1 == guppy_call ] ; then
    mkdir -p $ssddir/called
    for i in $samps ;
    do
	mkdir -p $ssddir/called/$i
	guppy_basecaller \
	    -i $ssddir/raw/$i \
	    -s $ssddir/called/$i \
	    --recursive \
	    --compress_fastq \
	    --flowcell FLO-MIN106 --kit SQK-LSK108 \
	    --device 'cuda:0'
    done
fi

if [ $1 == gather ] ; then
    mkdir -p $datadir/fastqs
    for i in $samps ;
    do
	mkdir -p $datadir/fastqs/$i
	cat $ssddir/called/$i/pass/*.fastq.gz > $datadir/fastqs/$i/$i.fq.gz
    done
fi

if [ $1 == flye ] ; then
    mkdir -p $datadir/flye
    for i in $samps ;
    do
	mkdir -p $datadir/flye/$i
	flye \
	    --nano-raw $datadir/fastqs/$i/$i.fq.gz \
	    -o $datadir/flye/$i \
	    -t 36 \
	    -g 100m \
	    --plasmids \
	    --meta
	mv $datadir/flye/$i/assembly.fasta $datadir/flye/$i/$i.assembly.fasta
    done
fi


if [ $1 == megalodon_contigs ] ; then
    mkdir -p $ssddir/megalodon/${i}_meta
    for i in $samps ;
    do
	asm=$datadir/flye/$i/$i.assembly.fasta
	echo $i
	megalodon \
	    $ssddir/raw/$i \
	    --overwrite \
	    --suppress-progress-bars \
	    --verbose-read-progress 0 \
	    --guppy-server-path "/usr/bin/guppy_basecall_server" \
	    --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	    --reference $asm \
	    --outputs per_read_mods \
	    --output-directory $ssddir/megalodon/${i}_meta \
	    --write-mods-text \
	    --devices "cuda:0" \
	    --processes 36
	rm $ssddir/megalodon/${i}_meta/per_read_modified_base_calls.db
    done
fi
	
dbdir=/mithril/Data/Nanopore/ref/kraken2

if [ $1 == kraken ] ; then
    mkdir -p $datadir/kraken
    for i in $samps ;
    do
	kraken2 \
	    --db $dbdir/standard_ont \
	    --threads 36 \
	    --classified-out $datadir/kraken/$i.class.txt \
	    --unclassified-out $datadir/kraken/$i.unclass.txt \
	    --output $datadir/kraken/$i.out.txt \
	    --report $datadir/kraken/$i.report.txt \
	    --use-names \
	    $datadir/fastqs/$i/$i.fq.gz
    done
fi


if [ $1 == top40 ] ; then
    for i in $samps ;
    do
	awk '$4 == "S" {print $0}' $datadir/kraken/$i.report.txt | \
	    sort -r -k2 -n | \
	    head -n 40 > $datadir/kraken/$i.report.top40.txt
    done
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
	rm $ssddir/megalodon/$i/per_read_modified_base_calls.db
    done
fi


if [ $1 == megaidx ] ; then
    for i in $samps ;
    do
	python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
		-i $ssddir/megalodon/$i/per_read_modified_base_calls.txt \
		-o $ssddir/megalodon/$i/per_read_modified_base_calls.txt.idx
    done
fi



if [ $1 == barcode ] ; then
    mkdir -p $ssddir/barcode
    for i in $samps ;
    do
	mkdir -p $ssddir/barcode/$i
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $ssddir/megalodon/$i/per_read_modified_base_calls.txt \
               -i $ssddir/megalodon/$i/per_read_modified_base_calls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt \
               -o $ssddir/barcode/$i/${i}_barcodes.txt \
               -t 36 ;} &> $ssddir/barcode/$i/${i}_time.txt
    done
fi


if [ $1 == copy ] ; then
    for i in $samps ;
    do
	cp -r $ssddir/barcode/$i $datadir/barcode/
	cp -r $ssddir/called/$i $datadir/called/
	cp -r $ssddir/megalodon/$i $datadir/megalodon/
    done
fi



if [ $1 == align ] ; then
    mkdir -p $datadir/align
    for i in $samps ;
    do
	mkdir -p $datadir/align/$i
	fq=$datadir/fastqs/$i/$i.fq.gz
	minimap2 -t 36 -x map-ont $ref $fq \
		 > $datadir/align/$i/$i.paf
    done
fi



if [ $1 == filter ] ; then
    for i in $samps ;
    do
	for j in 1 5 10 ;
	do
	    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
		   -a $datadir/align/$i/$i.paf \
		   -r $ref \
		   -o $datadir/barcode/$i/${i}_barcodes_filtered.$j.txt \
		   -m $datadir/barcode/$i/${i}_barcodes.txt \
		   -b ~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt \
		   -q 40 \
		   -l 5000 \
		   -n $j \
		   -v
	done
    done
fi

if [ $1 == remegaidx ] ; then
    for i in $samps ;
    do
	python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
		-i $datadir/megalodon/$i/per_read_modified_base_calls.txt \
		-o $datadir/megalodon/$i/per_read_modified_base_calls.txt.idx
    done
fi


if [ $1 == rebarcode ] ; then
    ##barocde on nas after fixing bugs 
    mkdir -p $datadir/barcode
    for i in $samps ;
    do
	mkdir -p $datadir/barcode/$i
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/$i/per_read_modified_base_calls.txt \
               -i $datadir/megalodon/$i/per_read_modified_base_calls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt \
               -o $datadir/barcode/$i/${i}_barcodes.txt \
               -t 36 ;} &> $datadir/barcode/$i/${i}_time.txt
    done
fi

if [ $1 == rebarcode_test ] ; then
    mkdir -p $datadir/barcode
    for i in MinION_BA_NAT ;
    do
	mkdir -p $datadir/barcode/$i
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/$i/per_read_modified_base_calls.txt \
               -i $datadir/megalodon/$i/per_read_modified_base_calls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt \
               -o $datadir/barcode/$i/${i}_barcodes_test.txt \
               -t 36 ;} &> $datadir/barcode/$i/${i}_time_test.txt
    done
fi

	
if [ $1 == barcode_common ] ; then
    mkdir -p $datadir/barcode
    for i in $samps ;
    do
	mkdir -p $datadir/barcode/$i
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/$i/per_read_modified_base_calls.txt \
               -i $datadir/megalodon/$i/per_read_modified_base_calls.txt.idx \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
               -o $datadir/barcode/$i/${i}_barcodes.txt \
               -t 12 ;} &> $datadir/barcode/$i/${i}_time15.txt
    done
fi


if [ $1 == megaidx_contig ] ; then
    for i in $samps ;
    do
	python ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	       -i $datadir/megalodon/${i}_meta/per_read_modified_base_calls.txt \
	       -o $datadir/megalodon/${i}_meta/per_read_modified_base_calls.txt.idx
    done
fi


if [ $1 == barcode_contig ] ; then
    mkdir -p $datadir/barcode
    for i in $samps ;
    do
	asm=$datadir/flye/$i/$i.assembly.fasta
	mkdir -p $datadir/barcode/${i}_meta
	{ time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
               -m $datadir/megalodon/${i}_meta/per_read_modified_base_calls.txt \
               -i $datadir/megalodon/${i}_meta/per_read_modified_base_calls.txt.idx \
               -r $asm \
               -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
               -o $datadir/barcode/${i}_meta/${i}_meta_barcodes15.txt \
               -t 12 ;} &> $datadir/barcode/${i}_meta/${i}_meta_time15.txt
    done
fi


if [ $1 == align ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/align/$i
	minimap2 -t 36 -x map-ont $ref $datadir/fastqs/$i/$i.fq.gz \
		 > $datadir/align/$i/$i.paf

	asm=$datadir/flye/$i/$i.assembly.fasta
	mkdir -p $datadir/align/${i}_meta
	minimap2 -t 36 -x map-ont $asm $datadir/fastqs/$i/$i.fq.gz \
		 > $datadir/align/${i}_meta/${i}_meta.paf
    done
fi

if [ $1 == motifcounts ] ; then
    for i in $samps ;
    do
	asm=$datadir/flye/$i/$i.assembly.fasta
	python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	       -a $datadir/align/$i/${i}.paf \
	       -r $ref \
	       -o $datadir/barcode/$i/${i}_motifcounts.txt \
	       -m $datadir/megalodon/$i/${i}_barcodes.txt \
	       -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
	       -q 40 \
	       -v
	python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	       -a $datadir/align/${i}_meta/${i}_meta.paf \
	       -r $asm \
	       -o $datadir/barcode/${i}_meta/${i}_meta_motifcounts.txt \
	       -m $datadirx/barcode/${i}_meta/${i}_meta_barcodes15.txt \
	       -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
	       -q 40 \
	       -v
    done
fi

blastdbdir=/atium/Data/ref/bacteria/blastdb
if [ $1 == blast_contigs ] ; then
    mkdir -p $datadir/blast_contigs

    i=MinION_JM3O_NAT
    blastn \
	-query $datadir/flye/$i/$i.assembly.fasta \
	-db $blastdbdir/all_bacteria_refs \
	-outfmt 7 \
	-num_threads 36 \
	-out $datadir/blast_contigs/flye_blast.tsv
fi
