#!/bin/bash

##This time, $1 is the pilon outdir. $2 is the assembly full path. $3 is the prefix. 

mkdir -p $1
mkdir -p $1/btidx
mkdir -p $1/btbam


##copy the raw assembly into the index dir
prefix=$3
cp $2 $1/btidx/$prefix.fasta
cat $1/*fq.gz > $1/all.fq.gz

for i in {1..15} ;
do
    ##build the index and align
    echo building btidx and aligning for round $i
    bowtie2-build -q $1/btidx/$prefix.fasta $1/btidx/$prefix
    bowtie2 -p 24 -x $1/btidx/$prefix -U $1/all.fq.gz > $1/btbam/$prefix.sam

    ##do the correction
    echo correcting for round $i
    racon -t 36 $1/all.fq.gz $1/btbam/$prefix.sam $1/btidx/$prefix.fasta > $1/$prefix.racon.$i.fasta

    ##newly corrected genome replaces the old genome in the index dir
    echo clearing old
    rm $1/btidx/*
    rm $1/btbam/*
    echo copying $i to empty index folder
    cp $1/$prefix.racon.$i.fasta $1/btidx/$prefix.fasta
done
    ##sed -i -e 's/_pilon//g' $1/pilon/btidx/$prefix.fasta

