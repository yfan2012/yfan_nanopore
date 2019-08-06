#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/jose_burkle

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    tarpath=/work-zfs/mschatz1/cpowgs/jose_burkle/190117_jose_burkle.tar.gz
    python3 ~/Code/timp_nanopore/aws/tar4kbin.py -i $tarpath --outdir $datadir/raw
fi


if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done

    bash ./bc_call.sh $datadir
fi

if [ $1 == fqs ] ; then
    mkdir -p $datadir/fastqs

    cat $datadir/called/*/workspace/pass/barcode01/*.fastq > $datadir/fastqs/jose1.fq
    cat $datadir/called/*/workspace/pass/barcode02/*.fastq > $datadir/fastqs/jose2.fq
    cat $datadir/called/*/workspace/pass/barcode03/*.fastq > $datadir/fastqs/jose3.fq
fi

if [ $1 == assemble ] ; then
    ##assemble with canu
    ##for i in jose1 jose2 jose3 ;
    for i in jose1 jose2 ;
    do
	mkdir -p $datadir/${i}_assembly
	python ~/Code/utils/fastq_long.py -i $datadir/fastqs/$i.fq -o $datadir/fastqs/${i}_3kb.fq -l 3000
	canu \
	    -p $i -d $datadir/${i}_assembly \
	    -gridOptions="--time=22:00:00 --partition=parallel" \
	    genomeSize=8.5m \
	    stopOnReadQuality=false \
	    -nanopore-raw $datadir/fastqs/${i}_3kb.fq
    done
fi

if [ $1 == trim ] ; then
    ##trim illumina reads 
    mkdir -p $datadir/illumina_trimmed
    for i in jose1 jose2 jose3 ;
    do
	java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 \
	     $datadir/illumina/${i}*R1*.gz $datadir/illumina/${i}*R2*.gz \
	     $datadir/illumina_trimmed/${i}_forward_paired.fq.gz $datadir/illumina_trimmed/${i}_forward_unpaired.fq.gz \
	     $datadir/illumina_trimmed/${i}_reverse_paired.fq.gz $datadir/illumina_trimmed/${i}_reverse_unpaired.fq.gz \
	     ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

	
if [ $1 == circlate ] ; then
    ##circlate
    ml samtools
    ml python
    ml gcc
    ##for i in jose1 jose2 jose3 ;
    for i in jose1 jose2 jose3 ;
    do
	circlator all --assembler canu --merge_min_id 85 --merge_breaklen 1000 $datadir/${i}_assembly/*contigs.fasta $datadir/${i}_assembly/*trimmedReads.fasta.gz $datadir/${i}_circlator
    done
fi

if [ $1 == pilon ] ; then
    for i in jose3 ;
    do
	mkdir -p $datadir/${i}_pilon
	cp $datadir/illumina_trimmed/$i* $datadir/${i}_pilon
	bash ./pilon.sh $datadir/${i}_pilon $datadir/${i}_circlator/06.fixstart.fasta $i canucirc
    done
fi

	
