#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans


if [ $1 == transposition ] ; then
    mumdir=$datadir/compare/st31_vs_st90853/mummer

    mkdir -p $datadir/compare
    mkdir -p $mumdir
    asm1=$datadir/st31/final/st31.final.fasta
    asm2=$datadir/st90853/final/st90853.final.fasta

    nucmer \
	-p $mumdir/st31_vs_st90853 \
	$asm1 \
	$asm2
    mummerplot --filter --fat --postscript \
	       -p $mumdir/st31_vs_st90853 \
	       $mumdir/st31_vs_st90853.delta
    mummerplot --filter --fat --png \
	       -p $mumdir/st31_vs_st90853 \
	       $mumdir/st31_vs_st90853.delta
    dnadiff \
	-p $mumdir/st31_vs_st90853 \
	$asm1 \
	$asm2
fi
    
if [ $1 == crossmap ] ; then
    ##cross map reads between samples to see breakpoints
    mkdir -p $datadir/compare/st31_vs_st90853/align

    st31fa=$datadir/st31/final/st31.final.fasta
    st31reads=$datadir/st31/reads/ont/st31_long.fastq.gz 
    st90853fa=$datadir/st90853/final/st90853.final.fasta
    st90853reads=$datadir/st90853/reads/ont/st90853_long.fastq.gz
    
    minimap2 -t 36 -ax map-ont $st31fa $st90853reads |
        samtools view -@ 36 -b |
        samtools sort -@ 36 -o $datadir/compare/st31_vs_st90853/align/st31.st90853reads.sorted.bam
    samtools index $datadir/compare/st31_vs_st90853/align/st31.st90853reads.sorted.bam

    minimap2 -t 36 -ax map-ont $st90853fa $st31reads |
        samtools view -@ 36 -b |
        samtools sort -@ 36 -o $datadir/compare/st31_vs_st90853/align/st90853.st31reads.sorted.bam
    samtools index $datadir/compare/st31_vs_st90853/align/st90853.st31reads.sorted.bam
fi
    
if [ $1 == last_mummer ] ; then
    ##compare everything to the perfect st90853 to see if any broken tigs can be joined

    ref=$datadir/st90853/final/st90853.final2.fasta
    
    for i in st31.final st5317.final2 ;
    do
	strain=`echo $i | cut -d "." -f 1`
	mumdir=$datadir/$strain/mummer_against_st90853
	mkdir -p $mumdir
	
	nucmer \
	    -p $mumdir/${strain}_against_st90853 \
	    $ref \
	    $datadir/$strain/final/$i.fasta
	
	mummerplot --filter --fat --postscript \
		   -p $mumdir/${strain}_against_st90853 \
		   $mumdir/${strain}_against_st90853.delta
	mummerplot --filter --fat --png \
		   -p $mumdir/${strain}_against_st90853 \
		   $mumdir/${strain}_against_st90853.delta
	dnadiff \
	    -p $mumdir/${strain}_against_st90853 \
	    $ref \
	    $datadir/$strain/final/$i.fasta
    done
fi

    
