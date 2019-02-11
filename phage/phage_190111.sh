#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/phage
srcdir=~/Code/utils/marcc

##just taking basecalled fastqs from gridion

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
    mkdir -p $datadir/190111_phage/align
    rm $datadir/190111_phage/align/*
    for i in $datadir/190111_phage/fastqs/*.fq ;
    do
	prefix=` basename $i .fq `
	minimap2 -a -x map-ont -t 36 $datadir/Mycobacteriophages-All.fasta $i | samtools view -b | samtools sort -o $datadir/190111_phage/align/$prefix.sorted.bam -T $datadir/190111_phage/align/$prefix.reads.tmp - 
	samtools index $datadir/190111_phage/align/$prefix.sorted.bam
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
    ##datadir=/kyber/Data/Nanopore/phage/align
    ml samtools
    for i in $datadir/190111_phage/align/*.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	samtools view -b -F 0x100 $i | samtools sort -o $datadir/190111_phage/align/$prefix.primary.sorted.bam
	samtools index $datadir/190111_phage/align/$prefix.primary.sorted.bam
    done
fi

	     
if [ $1 == count ] ; then
    ##outdir=~/Dropbox/Timplab_Data/phage/counts
    outdir=$datadir/190111_phage/counts
    ref=$datadir/Mycobacteriophages-All.fasta
    for i in $datadir/190111_phage/align/*.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	python ~/Code/yfan_nanopore/phage/genome_counts.py -i $i -o $outdir/$prefix.csv -r $ref
    done
fi

	
