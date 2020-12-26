#!/bin/bash

datadir=/kyber/Data/Nanopore/projects/nina
asmdir=$datadir/assemblies/pilon

if [ $1 == mummer ] ; then
    ##mum all the things together
    mumdir=$datadir/mummer

    for i in $asmdir/*fasta ;
    do
        prefix=`basename $i .fasta | cut -d . -f 1`
	mkdir -p $mumdir/$prefix

	for j in $asmdir/*fasta ;
	do
	    qfix=`basename $j .fasta | cut -d . -f 1`
	    nucmer -t 36 -p $mumdir/$prefix/$qfix $i $j
	    mummerplot --png --fat --filter -p $mumdir/$prefix/$qfix.layout $mumdir/$prefix/$qfix.delta -R $i -Q $j
	    dnadiff -p $mumdir/$prefix/$qfix -d $mumdir/$prefix/$qfix.delta
	done
    done
fi


if [ $1 == minimap ] ; then
    alndir=$datadir/align

    for i in $asmdir/*fasta ;
    do
	asmprefix=`basename $i .fasta | cut -d . -f 1`
	prefix=`basename $i .fasta | cut -d _ -f 1`
	echo $prefix

	minimap2 -t 36 -ax map-ont $i $datadir/reads/$prefix/*over3kb.fastq | samtools view -@ 36 -b | samtools sort -@ 36 -o $alndir/$asmprefix.sorted.bam
	samtools index $alndir/$asmprefix.sorted.bam
    done
fi

	
if [ $1 == filtprimary ] ; then
    alndir=$datadir/align

    for i in $alndir/*sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	samtools view -@ 36 -b -F 2304 $i | samtools sort -@ 36 -o $alndir/$prefix.primary.sorted.bam
	samtools index $alndir/$prefix.primary.sorted.bam
    done
fi

if [ $1 == idxfa ] ; then
    for i in $asmdir/*fasta ;
    do
	prefix=`basename $i .sorted.bam`
	samtools faidx $i
    done
fi


if [ $1 == align_illumina ] ; then
    alndir=$datadir/align_illumina
    mkdir -p $alndir
    for i in $asmdir/*fasta ;
    do
	asmprefix=`basename $i .fasta | cut -d . -f 1`
	prefix=`basename $i .fasta | cut -d _ -f 1`
	echo $prefix

	samtools faidx $i
	
	bwa index $i
	cat $datadir/reads/$prefix/*.fq.gz > $datadir/reads/$prefix/all.fq.gz
	bwa mem -t 36 $i $datadir/reads/$prefix/all.fq.gz | samtools view -@ 36 -b | samtools sort -@ 36 -o $alndir/$asmprefix.sorted.bam
	samtools index $alndir/$asmprefix.sorted.bam
    done
fi
