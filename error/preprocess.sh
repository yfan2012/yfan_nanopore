#!/bin/bash

##script for everything up to assembly
datadir=/uru/Data/Nanopore/projects/nivar

if [ $1 == r10_unpack ] ; then
    prefix=r10
    
    tarball=/kyber/Data/Nanopore/oxford/190706_nivar_r10/190706_nivar_r10.tar.gz

    rawdir=$datadir/$prefix/raw
    mkdir -p $rawdir
    
    tar -xzf $tarball -C $datadir/$prefix/raw
fi


if [ $1 == r9_call ] ; then
    ##raw=/uru/Data/Nanopore/projects/nivar/r9/raw
    ##mounted /uru/Data/Nanopore/projects/nivar/r9 to /media
    raw=/media/r9/raw
    guppy_basecaller -i $raw/run1 -s $raw/run1_called --flowcell FLO-MIN106 --kit SQK-LSK109 -x "cuda:0"
    guppy_basecaller -i $raw/run2 -s $raw/run2_called --flowcell FLO-MIN106 --kit SQK-LSK109 -x "cuda:0"
fi


if [ $1 == r10_call ] ; then
    raw=/media/r10/raw/190706_nivar_r10
    guppy_basecaller -i $raw/fast5 -s $raw/called --flowcell FLO-MIN110 --kit SQK-LSK109 -x "cuda:0"
fi


if [ $1 == dRNA_call ] ; then
    raw=/media/dRNA
    guppy_basecaller -i $raw/dRNA1 -s $raw/dRNA1_called --flowcell FLO-MIN106 --kit SQK-RNA002 -x "cuda:0"
    guppy_basecaller -i $raw/dRNA2 -s $raw/dRNA2_called --flowcell FLO-MIN106 --kit SQK-RNA002 -x "cuda:0"
fi

##rearranged some of the called/raw data
if [ $1 == gather_fqs ] ; then
    dRNA=$datadir/dRNA/dRNA.fq
    r10=$datadir/r10/r10.fq
    r9=$datadir/r9/r9.fq

    cat $datadir/dRNA/dRNA1_called/*fastq $datadir/dRNA/dRNA2_called/*fastq > $dRNA
    cat $datadir/r9/run1_called/*fastq $datadir/r9/run2_called/*fastq > $r9
    cat $datadir/r10/called/*fastq > $r10

    python ~/Code/utils/fastq_long.py -i $r9 -o $datadir/r9/r9_3kb.fq -l 3000
    python ~/Code/utils/fastq_long.py -i $r10 -o $datadir/r10/r10_3kb.fq -l 3000
fi

if [ $1 == assemble ] ; then
    mkdir -p $datadir/assemble

    for i in r9 r10 ;
    do
	mkdir -p $datadir/assemble/${i}_assembly
	canu \
	    -p nivar_${i} -d $datadir/assemble/${i}_assembly \
	    genomeSize=12m \
	    -nanopore-raw $datadir/$i/${i}_3kb.fq
    done
fi

if [ $1 == trim ] ; then
    for i in cDNA gDNA ;
    do
	cat $datadir/illumina/$i/*_R1_*fastq.gz > $datadir/illumina/$i/nivar_${i}_R1.fastq.gz
	cat $datadir/illumina/$i/*_R2_*fastq.gz > $datadir/illumina/$i/nivar_${i}_R2.fastq.gz

	mkdir -p $datadir/illumina/${i}_trimmed
	
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 54 -phred33 \
	     $datadir/illumina/$i/nivar_${i}_R1.fastq.gz $datadir/illumina/$i/nivar_${i}_R2.fastq.gz \
	     $datadir/illumina/${i}_trimmed/nivar_${i}_fwd_paired.fq.gz $datadir/illumina/${i}_trimmed/nivar_${i}_fwd_unpaired.fq.gz \
	     $datadir/illumina/${i}_trimmed/nivar_${i}_rev_paired.fq.gz $datadir/illumina/${i}_trimmed/nivar_${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:idt_ud_indexes.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == freebayes ] ; then
    mkdir -p $datadir/freebayes
    for i in r9 r10 ;
    do
	mkdir -p $datadir/freebayes/${i}_freebayes
	cp $datadir/illumina/gDNA_trimmed/nivar_gDNA_*_paired.fq.gz $datadir/freebayes/${i}_freebayes/
	bash ~/Code/yfan_nanopore/nivar/scripts/freebayes_bwa.sh $datadir/freebayes/${i}_freebayes $datadir/assemble/${i}_assembly/nivar_$i.contigs.fasta nivar_$i
    done
fi

if [ $1 == pilon ] ; then
    mkdir -p $datadir/pilon
    for i in r9 r10 ;
    do
	mkdir -p $datadir/pilon/${i}_pilon
        cp $datadir/illumina/gDNA_trimmed/nivar_gDNA_*_paired.fq.gz $datadir/pilon/${i}_pilon
	bash ~/Code/yfan_nanopore/nivar/scripts/pilon_bwa.sh $datadir/pilon/${i}_pilon $datadir/assemble/${i}_assembly/nivar_$i.contigs.fasta nivar_$i
    done
fi

if [ $1 == clean_pilon ] ; then
    for i in r9 r10 ;
    do
	sed -i -e 's/_pilon//g' $datadir/pilon/${i}_pilon/nivar*15*
	sed -i -e 's/>/>pilon_15_/g' $datadir/pilon/${i}_pilon/nivar*15*
    done
fi


if [ $1 == racon ] ; then
    mkdir -p $datadir/racon
    for i in r9 r10 ;
    do
	mkdir -p $datadir/racon/${i}_racon
	cp $datadir/illumina/gDNA_trimmed/nivar_gDNA_*_paired.fq.gz $datadir/racon/${i}_racon
	bash ~/Code/yfan_nanopore/nivar/scripts/racon_bwa.sh $datadir/racon/${i}_racon $datadir/assemble/${i}_assembly/nivar_$i.contigs.fasta nivar_$i
    done
fi


if [ $1 == mummer ] ; then
    ref=$datadir/reference/candida_nivariensis.fa
    tr '[:lower:]' '[:upper:]' < $datadir/reference/candida_nivariensis_original.fa > $ref
    
    mkdir -p $datadir/mummer

    for i in r9 r10
    do
	raw=$datadir/assemble/${i}_assembly/nivar_${i}.contigs.fasta

	
	cp $ref ~/tmp/
	cp $raw ~/tmp/
	
	nucmer -p ~/tmp/nivar_${i}_raw_ref ~/tmp/candida_nivariensis.fa ~/tmp/nivar_${i}.contigs.fasta
	mummerplot --filter --fat --png -p ~/tmp/nivar_${i}_raw_ref ~/tmp/nivar_${i}_raw_ref.delta
	dnadiff -p ~/tmp/nivar_${i}_raw_ref ~/tmp/candida_nivariensis.fa ~/tmp/nivar_${i}.contigs.fasta


	mkdir -p $datadir/mummer/${i}_raw_ref
	cp ~/tmp/nivar* $datadir/mummer/${i}_raw_ref/
	
	rm ~/tmp/*

	
	for corr in pilon racon freebayes ;
	do
	    
	    ##mummer against raw
	    mkdir -p $datadir/mummer/${i}_${corr}_raw
	    
	    cp $raw ~/tmp/
	    cp $datadir/$corr/${i}_$corr/*15*.fasta ~/tmp/$corr.15.raw.fasta
	    tr '[:lower:]' '[:upper:]' < ~/tmp/$corr.15.raw.fasta > ~/tmp/$corr.15.fasta
	    sed -i -e "s/>/>${i}_${corr}_15_/g" ~/tmp/$corr.15.fasta
	    
	    
	    nucmer -p ~/tmp/nivar_${i}_${corr}_raw ~/tmp/$corr.15.fasta ~/tmp/nivar_${i}.contigs.fasta 
	    mummerplot --filter --fat --png -p ~/tmp/nivar_${i}_${corr}_raw ~/tmp/nivar_${i}_${corr}_raw.delta
	    dnadiff -p ~/tmp/nivar_${i}_${corr}_raw ~/tmp/$corr.15.fasta ~/tmp/nivar_${i}.contigs.fasta 
	    
	    cp ~/tmp/* $datadir/mummer/${i}_${corr}_raw/

	    rm ~/tmp/*


	    #mummer against reference
	    mkdir -p $datadir/mummer/${i}_${corr}_ref
	    
	    cp $ref ~/tmp/
	    cp $datadir/$corr/${i}_$corr/*15*.fasta ~/tmp/$corr.15.raw.fasta
	    tr '[:lower:]' '[:upper:]' < ~/tmp/$corr.15.raw.fasta > ~/tmp/$corr.15.fasta
	    sed -i -e "s/>/>${i}_${corr}_15_/g" ~/tmp/$corr.15.fasta
	    
	    nucmer -p ~/tmp/nivar_${i}_${corr}_ref ~/tmp/candida_nivariensis.fa ~/tmp/$corr.15.fasta
	    mummerplot --filter --fat --png -p ~/tmp/nivar_${i}_${corr}_ref ~/tmp/nivar_${i}_${corr}_ref.delta
	    dnadiff -p ~/tmp/nivar_${i}_${corr}_ref ~/tmp/candida_nivariensis.fa ~/tmp/$corr.15.fasta

	    cp ~/tmp/* $datadir/mummer/${i}_${corr}_ref/
	    
	    rm ~/tmp/*
	    ##mummer fb against the reference
	    
	done
    done
fi


	
    
