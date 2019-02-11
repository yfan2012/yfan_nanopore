#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/moth

if [ $1 == assemble ] ; then
   wtdbg2 -t 32 -i $datadir/181005_moth_cat.fastq -fo $datadir/wtdbg2
fi

if [ $1 == cons ] ; then
   wtpoa-cns -t 32 -i $datadir/wtdbg2.ctg.lay.gz -fo $datadir/wtdbg2.ctg.lay.fasta
fi

if [ $1 == canu ] ; then
    canu \
	-p moth -d $datadir/canu \
	-gridOptions="--time=72:00:00 --partition=parallel" \
	genomeSize=500m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/181005_moth_cat.fastq
fi

if [ $1 == alignpilon ] ; then
    ml samtools
    asm=$datadir/wtdbg2_pilon/moth.wtdbg2.pilon.13.fasta
    npreads=$datadir/181005_moth_cat.fastq
    
    minimap2 -a -x map-ont -t 36 $asm $npreads | samtools view -b | samtools sort -o $datadir/align/moth.wtdbg2.pilon.13.npreads.sorted.bam -T $datadir/align/reads.tmp -
    samtools index $datadir/align/moth.wtdbg2.pilon.13.npreads.sorted.bam

    mkdir $datadir/align/btidx
    bowtie2-build -q $asm $datadir/align/btidx/moth.wtdbg2.pilon.13
    bowtie2 -p 36 -x $datadir/align/btidx/moth.wtdbg2.pilon.13 -1 $datadir/*R1_001.fastq.gz -2 $datadir/*R2_001.fastq.gz | samtools view -bS - | samtools sort -o $datadir/align/moth.wtdbg2.pilon.13.illreads.sorted.bam
    samtools index $datadir/align/moth.wtdbg2.pilon.13.illreads.sorted.bam
fi
    
