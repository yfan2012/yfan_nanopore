#!/bin/bash

phagedir=/scratch/groups/mschatz1/cpowgs/phage
datadir=$phagedir/190412_phage
														      

if [ $1 == align ] ; then
    ml samtools

    mkdir -p $datadir/align

    for i in $datadir/fastqs/*fastq ;
    do
	prefix=`basename $i .fastq`
	echo $prefix
	minimap2 -a -x map-ont -t 36 $phagedir/Mycobacteriophages-All.fasta $i | samtools view -b | samtools sort -o $datadir/align/$prefix.sorted.bam -T $datadir/align/$prefix.reads.tmp
	samtools index $datadir/align/$prefix.sorted.bam
    done
fi

if [ $1 == counts ] ; then
    outdir=$datadir/counts
    mkdir -p $outdir

    ref=$datadir/Mycobacteriophages-All.fasta
    for i in $datadir/align/*.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	python ~/Code/yfan_nanopore/phage/genome_counts.py -i $i -o $outdir/$prefix.csv -r $ref
    done
fi
