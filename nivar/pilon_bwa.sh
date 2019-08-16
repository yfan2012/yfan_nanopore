#!/bin/bash

##This time, $1 is the pilon outdir. $2 is the assembly full path. $3 is the prefix. 

mkdir -p $1
mkdir -p $1/index
mkdir -p $1/bam


##copy the raw assembly into the index dir
prefix=$3
cp $2 $1/index/$prefix.fasta


for i in {1..15} ;
do
    ##build the index and align
    echo building index and aligning for round $i
    bwa index $1/index/$prefix.fasta
    bwa mem -t 36 $1/index/$prefix.fasta $1/*forward_paired.fq.gz $1/*reverse_paired.fq.gz | samtools view -bS - | samtools sort -o $1/bam/$prefix.sorted.bam
    samtools index $1/bam/$prefix.sorted.bam

    ##do the correction
    echo correcting for round $i
    java -Xmx100G -jar ~/software/pilon-1.23.jar --threads 12 --changes --tracks --genome $1/index/$prefix.fasta --frags $1/bam/$prefix.sorted.bam --outdir $1 --output $prefix.pilon_bwa.$i

    ##newly corrected genome replaces the old genome in the index dir
    echo clearing old
    rm $1/index/*
    rm $1/bam/*
    echo copying $i to empty index folder
    cp $1/$prefix.pilon_bwa.$i.fasta $1/index/$prefix.fasta
done


