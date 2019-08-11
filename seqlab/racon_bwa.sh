#!/bin/bash

##This time, $1 is the pilon outdir. $2 is the assembly full path. $3 is the prefix. 

mkdir -p $1
mkdir -p $1/index

##copy the raw assembly into the index dir
prefix=$3
cp $2 $1/index/$prefix.fasta
cat $1/*fq.gz > $1/all.fq.gz

for i in {1..15} ;
do
    ##build the index and align
    echo building btidx and aligning for round $i
    bwa index $1/index/$prefix.fasta
    bwa mem -t 36 $1/index/$prefix.fasta $1/all.fq.gz > $1/$prefix.sam

    ##do the correction
    echo correcting for round $i
    racon -t 36 $1/all.fq.gz $1/$prefix.sam $1/index/$prefix.fasta > $1/$prefix.racon.$i.fasta

    ##newly corrected genome replaces the old genome in the index dir
    echo clearing old
    rm $1/index/*
    echo copying $i to empty index folder
    cp $1/$prefix.racon.$i.fasta $1/index/$prefix.fasta
done
    ##sed -i -e 's/_pilon//g' $1/pilon/btidx/$prefix.fasta

