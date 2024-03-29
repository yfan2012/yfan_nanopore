#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/phage
srcdir=~/Code/utils/marcc

if [ $1 == untar ] ; then
    mkdir -p $datadir/batch_logs
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_phage $srcdir/untar.scr $datadir/180809_bacteriophage_barcoded_15_20.raw.tgz $datadir
fi
    
if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done
    sbatch --array=0-149 --job-name=call_phage --output=$datadir/call_logs/phage.%A_%a.out $srcdir/bc_call.scr $datadir
fi

if [ $1 == fastqs ] ; then
    mkdir -p $datadir/fastqs
    cat $datadir/called/*/workspace/pass/barcode01/*.fastq > $datadir/fastqs/samp15at24hrs.fq
    cat $datadir/called/*/workspace/pass/barcode02/*.fastq > $datadir/fastqs/samp15at48hrs.fq
    cat $datadir/called/*/workspace/pass/barcode03/*.fastq > $datadir/fastqs/samp20at24hrs.fq
    cat $datadir/called/*/workspace/pass/barcode04/*.fastq > $datadir/fastqs/samp20at48hrs.fq
fi


if [ $1 == vanilla_centrifuge ] ; then
    ml gcc
    dbdir=~/scratch/centrifuge_db
    mkdir -p $datadir/classification
    for i in $datadir/fastqs/*.fq ;
    do
	prefix=`basename $i .fq`
	##~/scratch/centrifuge/centrifuge -p 36 -x $dbdir/abvm -U $i -S $datadir/classification/${prefix}_classification.txt --report-file $datadir/classification/${prefix}_report.tsv
	~/scratch/centrifuge/centrifuge-kreport -x $dbdir/abvm $datadir/classification/${prefix}_classification.txt > $datadir/classification/${prefix}_kreport.txt
    done
fi


if [ $1 == align ] ; then
    ml samtools
    mkdir -p $datadir/180809_phage/align
    rm $datadir/180809_phage/align/*
    for i in $datadir/180809_phage/fastqs/*.fq ;
    do
	prefix=` basename $i .fq `
	minimap2 -a -x map-ont -t 36 $datadir/Mycobacteriophages-All.fasta $i | samtools view -b | samtools sort -o $datadir/180809_phage/align/$prefix.sorted.bam -T $datadir/180809_phage/align/$prefix.reads.tmp - 
	samtools index $datadir/180809_phage/align/$prefix.sorted.bam
    done
fi

    

if [ $1 == unique ] ; then
    datadir=/kyber/Data/Nanopore/phage/align
    for i in $datadir/*hrs.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	python ~/Code/yfan_nanopore/phage/exactly_one.py -i $i -o $datadir/$prefix.unique.bam
	samtools sort $datadir/$prefix.unique.bam -o $datadir/$prefix.unique.sorted.bam
	samtools index $datadir/$prefix.unique.sorted.bam
    done
fi

if [ $1 == primary ] ; then
    datadir=/kyber/Data/Nanopore/phage/align
    for i in $datadir/*hrs.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	samtools view -b -F 0x100 $i | samtools sort -o $datadir/$prefix.primary.sorted.bam
	samtools index $datadir/$prefix.primary.sorted.bam
    done
fi

	     
if [ $1 == count ] ; then
    ##outdir=~/Dropbox/Timplab_Data/phage/counts
    outdir=$datadir/180809_phage/counts
    ref=$datadir/Mycobacteriophages-All.fasta
    for i in $datadir/180809_phage/align/*.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	python ~/Code/yfan_nanopore/phage/genome_counts.py -i $i -o $outdir/$prefix.csv -r $ref
    done
fi

	
