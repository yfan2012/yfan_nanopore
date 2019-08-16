#!/bin/bash

datadir=~/work/seqlab/fungus

if [ $1 == long ] ; then
    python ~/Code/utils/fastq_long.py -i $datadir/candida_nivariensis.fastq -o $datadir/candida_nivariensis_3k.fastq -l 3000
fi

if [ $1 == canu ] ; then
    mkdir -p $datadir/canu
    canu \
	-p candida_nivariensis -d $datadir/canu/ \
	-gridOptions="--time=22:00:00 --partition=parallel" \
	genomeSize=15m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/fastqs/candida_nivariensis_3k.fastq
fi


if [ $1 == wtdbg2 ] ; then
    mkdir -p $datadir/wtdbg2

    wtdbg2 -t 32 -i $datadir/fastqs/candida_nivariensis_3k.fastq -fo $datadir/wtdbg2/candida_nivariensis_wtdbg2
    wtpoa-cns -t 32 -i $datadir/wtdbg2/candida_nivariensis_wtdbg2.ctg.lay.gz -fo $datadir/wtdbg2/candida_nivariensis.wtdbg2.contigs.fasta
fi


if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed

    for i in cDNA gDNA ;
    do
	java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 \
	     $datadir/${i}_illumina/CANI_${i}_R1.fastq.gz $datadir/${i}_illumina/CANI_${i}_R2.fastq.gz \
	     $datadir/trimmed/CANI_${i}_forward_paired.fq.gz $datadir/trimmed/CANI_${i}_forward_unpaired.fq.gz \
	     $datadir/trimmed/CANI_${i}_reverse_paired.fq.gz $datadir/trimmed/CANI_${i}_reverse_unpaired.fq.gz \
	     ILLUMINACLIP:idt_ud_indexes.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == pilon ] ; then
    mkdir -p $datadir/canu_pilon
    mkdir -p $datadir/wtdbg2_pilon

    ##cp $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz $datadir/canu_pilon
    ##cp $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz $datadir/canu_pilon
    cp $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz $datadir/wtdbg2_pilon
    cp $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz $datadir/wtdbg2_pilon
    
    
    ##bash ./pilon2.scr $datadir/canu_pilon $datadir/canu/candida_nivariensis.contigs.fasta candida_nivariensis canu
    bash ./pilon2.scr $datadir/wtdbg2_pilon $datadir/wtdbg2/candida_nivariensis.wtdbg2.contigs.fasta candida_nivariensis wtdbg2
fi

if [ $1 == spades ] ; then
    mkdir -p $datadir/spades

    spades.py -1 $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz -2 $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz -t 36 -m 300 -o $datadir/spades
fi
