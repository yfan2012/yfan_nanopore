#!/bin/bash

##This time, $1 is the pilon outdir. $2 is the assembly full path. $3 is the prefix. $4 is the assembler

ml samtools

mkdir -p $1
mkdir -p $1/btidx
mkdir -p $1/btbam


##copy the raw assembly into the index dir
prefix=$3_$4
cp $2 $1/btidx/$prefix.fasta


for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ;
do
    ##build the index and align
    echo building btidx and aligning for round $i
    bowtie2-build -q $1/btidx/$prefix.fasta $1/btidx/$prefix
    bowtie2 -p 24 -x $1/btidx/$prefix -1 $1/$3*forward_paired.fq.gz -2 $1/$3*reverse_paired.fq.gz | samtools view -bS - | samtools sort -o $1/btbam/$prefix.sorted.bam
    samtools index $1/btbam/$prefix.sorted.bam

    ##do the correction
    echo correcting for round $i
    java -Xmx100G -jar ~/software/pilon/pilon-1.22.jar --threads 12 --changes --tracks --genome $1/btidx/$prefix.fasta --frags $1/btbam/$prefix.sorted.bam --outdir $1 --output $prefix.pilon.$i

    ##newly corrected genome replaces the old genome in the index dir
    echo clearing old
    rm $1/btidx/*
    rm $1/btbam/*
    echo copying $i to empty index folder
    cp $1/$prefix.pilon.$i.fasta $1/btidx/$prefix.fasta
done
    ##sed -i -e 's/_pilon//g' $1/pilon/btidx/$prefix.fasta

