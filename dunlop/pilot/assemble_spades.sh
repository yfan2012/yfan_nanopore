#!/bin/bash

if [ $1 == wha??? ] ; then
    ml python/2.7
    
    ##set up some directory names
    datadir=~/work/nanopore
    fastqdir=$datadir/fastqs
    
    mkdir -p $datadir/assemble
    cp ~/software/Trimmomatic-0.38/adapters/NexteraPE-PE.fa ./
    
    
    ##do the following commands for all the samples in the directory
    for i in $fastqdir/*R1_001.fastq.gz;
    do
	##get the sample name from the file of reads
	prefix=`echo ${i#$fastqdir/} | cut -d '_' -f 1,2`
	echo $prefix
	
	mkdir -p $datadir/assemble/$prefix
	
	##assemble using spades
	
	java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 $i $fastqdir/${prefix}_L001_R2_001.fastq.gz $datadir/trimmed/${prefix}_forward_paired.fq.gz $datadir/trimmed/${prefix}_forward_unpaired.fq.gz $datadir/trimmed/${prefix}_reverse_paired.fq.gz $datadir/trimmed/${prefix}_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	
	spades.py -1 $datadir/trimmed/${prefix}_forward_paired.fq.gz -2 $datadir/trimmed/${prefix}_reverse_paired.fq.gz -t 36 -m 300 -o $datadir/assemble/$prefix
    done
    
fi


if [ $1 == old_data ] ; then 
    ml python/2.7.10
    
    ##set up some directory names
    datadir=~/work/ngs/180628_sears_fmt2
    fastqdir=$datadir/fastqs
    
    mkdir -p $datadir/assemble
    cp ~/software/Trimmomatic-0.38/adapters/NexteraPE-PE.fa ./
    
    
    ##do the following commands for all the samples in the directory
    for i in $fastqdir/*R1_001.fastq.gz;
    do
	##get the sample name from the file of reads
	prefix=`echo ${i#$fastqdir/} | cut -d '_' -f 1,2`
	echo $prefix
	
	mkdir -p $datadir/assemble/$prefix
	
	##assemble using spades
	
	java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 $i $fastqdir/${prefix}_L001_R2_001.fastq.gz $datadir/trimmed/${prefix}_forward_paired.fq.gz $datadir/trimmed/${prefix}_forward_unpaired.fq.gz $datadir/trimmed/${prefix}_reverse_paired.fq.gz $datadir/trimmed/${prefix}_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	
	spades.py -1 $datadir/trimmed/${prefix}_forward_paired.fq.gz -2 $datadir/trimmed/${prefix}_reverse_paired.fq.gz -t 36 -m 300 -o $datadir/assemble/$prefix
    done
fi


if [ $1 == mgrb ] ; then
    ml python/2.7

    ##compare illumina assembly vs nanopore for dunlop grant
    datadir=/scratch/groups/mschatz1/cpowgs/r21/KLPN_70
    fastqdir=$datadir/pilon
    prefix=KLPN_70
    
    mkdir -p $datadir/spades
    mkdir -p $datadir/spades/trimmed
    cp ~/software/Trimmomatic-0.38/adapters/NexteraPE-PE.fa ./

    java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 $fastqdir/*R1* $fastqdir/*_L001_R2_001.fastq.gz $datadir/spades/trimmed/${prefix}_forward_paired.fq.gz $datadir/spades/trimmed/${prefix}_forward_unpaired.fq.gz $datadir/spades/trimmed/${prefix}_reverse_paired.fq.gz $datadir/spades/trimmed/${prefix}_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    spades.py -1 $datadir/spades/trimmed/${prefix}_forward_paired.fq.gz -2 $datadir/spades/trimmed/${prefix}_reverse_paired.fq.gz -t 48 -m 300 -o $datadir/spades/$prefix
fi
