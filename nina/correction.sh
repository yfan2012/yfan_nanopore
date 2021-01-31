#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans

if [ $1 == align ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/racon
	fq=$datadir/$i/reads/ont/${i}_long.fastq.gz
	for asm in canu flye ;
	do
	    minimap2 -t 54 -ax map-ont $datadir/$i/asm/$asm/$i.c*.fasta $fq |
		samtools view -@ 54 -b |
		samtools sort -@ 54 -o $datadir/$i/racon/$i.$asm.sorted.bam
	    samtools index $datadir/$i/racon/$i.$asm.sorted.bam
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
	    bam=$datadir/$i/racon/$i.$asm.sorted.bam
	    samtools view -@ 54 $bam > $datadir/$i/racon/$i.$asm.sorted.sam
	    
	    ##settings recommended by medaka page
	    racon -m 8 -x -6 -g -8 -w 500 -t 54\
		  $fq \
		  $datadir/$i/racon/$i.$asm.sorted.sam \
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
    ##for i in st31 st5317 st90853 ;
    for i in st5317;
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
	    

if [ $1 == choose ] ; then
    ##manually inspecting vcfs to pick which ones have converged
    mkdir -p $datadir/st5317/genomes
    cp $datadir/st5317/freebayes/canu/st5317.canu_fb12.fasta $datadir/st5317/genomes/st5317.canu.fasta
    cp $datadir/st5317/freebayes/flye/st5317.flye_fb10.fasta $datadir/st5317/genomes/st5317.flye.fasta

    mkdir -p $datadir/st31/genomes
    cp $datadir/st31/freebayes/canu/st31.canu_fb3.fasta $datadir/st31/genomes/st31.canu.fasta
    cp $datadir/st31/freebayes/flye/st31.flye_fb4.fasta $datadir/st31/genomes/st31.flye.fasta

    mkdir -p $datadir/st90853/genomes
    cp $datadir/st90853/freebayes/canu/st90853.canu_fb4.fasta $datadir/st90853/genomes/st90853.canu.fasta
    cp $datadir/st90853/freebayes/flye/st90853.flye_fb5.fasta $datadir/st90853/genomes/st90853.flye.fasta
fi


if [ $1 == ragtag ] ; then

    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/ragtag/cf
	ragtag.py scaffold \
		  -w \
		  -u \
		  -o $datadir/$i/ragtag/cf \
		  $datadir/$i/genomes/$i.canu.fasta \
		  $datadir/$i/genomes/$i.flye.fasta

	mv $datadir/$i/ragtag/cf/ragtag.scaffolds.fasta $datadir/$i/ragtag/cf/$i.ragtag_cf.scaffolds.fasta

	mkdir -p $datadir/$i/ragtag/fc
	ragtag.py scaffold \
		  -w \
		  -u \
		  -o $datadir/$i/ragtag/fc \
		  $datadir/$i/genomes/$i.flye.fasta \
		  $datadir/$i/genomes/$i.canu.fasta


	mv $datadir/$i/ragtag/fc/ragtag.scaffolds.fasta $datadir/$i/ragtag/fc/$i.ragtag_fc.scaffolds.fasta

    done
fi
	    

    
if [ $1 == add_scaffolds ] ; then
    for i in st31 st5317 st90853 ;
    do
	cp $datadir/$i/ragtag/fc/$i.ragtag_fc.scaffolds.fasta $datadir/$i/genomes/$i.ragtag_fc.fasta
	cp $datadir/$i/ragtag/cf/$i.ragtag_cf.scaffolds.fasta $datadir/$i/genomes/$i.ragtag_cf.fasta
    done
fi
