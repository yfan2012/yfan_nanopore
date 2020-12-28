#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans

if [ $1 == align ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/align
	fq=$datadir/$i/reads/ont/${i}_long.fastq.gz
	for asm in canu flye ;
	do
	    minimap2 -t 54 -ax map-ont $datadir/$i/asm/$asm/$i.c*.fasta $fq |
		samtools view -@ 54 -b |
		samtools sort -@ 54 -o $datadir/$i/align/$i.$asm.sorted.bam
	    samtools index $datadir/$i/align/$i.$asm.sorted.bam
	done
    done
fi

if [ $1 == racon ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/racon
	fq=$datadir/$i/reads/ont/${i}_long.fastq.gz

	for asm in canu flye ;
	do
	    bam=$datadir/$i/align/$i.$asm.sorted.bam
	    samtools view -@ 54 $bam > $datadir/$i/align/$i.$asm.sorted.sam
	    
	    ##settings recommended by medaka page
	    racon -m 8 -x -6 -g -8 -w 500 -t 54\
		  $fq \
		  $datadir/$i/align/$i.$asm.sorted.sam \
		  $datadir/$i/asm/$asm/$i.c*.fasta > $datadir/$i/racon/$i.$asm.racon.contigs.fasta
	done
    done
fi
	      
	    

if [ $1 == medaka ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/medaka
	for asm in canu flye ;
	do

	    mkdir -p $datadir/$i/medaka/$asm
	    
	    medaka_consensus \
		-i $datadir/$i/reads/ont/${i}_long.fastq.gz \
		-d $datadir/$i/racon/$i.$asm.racon.contigs.fasta \
		-o $datadir/$i/medaka/$asm \
		-t 54 \
		-m r941_min_high_g360
	done
    done
fi

if [ $1 == freebayes ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/freebayes
	for asm in canu flye ;
	do

	    mkdir -p $datadir/$i/freebayes/$asm
	    cp $datadir/$i/medaka/$asm/consensus.fasta $datadir/$i/medaka/$asm/$i.$asm.medaka.fasta
	    
	    bash ~/Code/yfan_nanopore/nina/freebayes.sh \
		 $datadir/$i/freebayes/$asm \
		 $datadir/$i/medaka/$asm/$i.$asm.medaka.fasta \
		 $datadir/$i/reads/illumina \
		 $i.$asm
	done
    done
fi
	    
