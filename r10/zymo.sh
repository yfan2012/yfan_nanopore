#!/bin/bash

datadir=~/work/r10_zymo

if [ $1 == assembly ] ; then
    mkdir $datadir/canu_r10
    mkdir $datadir/canu_r9
    
    canu \
	-p zymo_r10 -d $datadir/canu_r10/ \
	-gridOptions="--time=22:00:00 --partition=parallel" \
	genomeSize=15m \
	corMinCoverage=0 \
	corOutCoverage=all \
	corMhapSensitivity=high \
	corMaxEvidenceCoverageLocal=10 \
	corMaxEvidenceCoverageGlobal=10 \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/fastqs/zymo_r10.fq.gz
    canu \
	-p zymo_r9 -d $datadir/canu_r9/ \
	-gridOptions="--time=22:00:00 --partition=parallel" \
	genomeSize=15m \
	corMinCoverage=0 \
	corOutCoverage=all \
	corMhapSensitivity=high \
	corMaxEvidenceCoverageLocal=10 \
	corMaxEvidenceCoverageGlobal=10 \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/fastqs/zymo_r9.fq.gz
fi

if [ $1 == ecoli_align ] ; then
    ml samtools

    ref=$datadir/Reference/ecoli_k12.fa
    
    mkdir -p $datadir/align_r10
    mkdir -p $datadir/align_r9

    minimap2 -a -x map-ont -t 36 $ref $datadir/fastqs/zymo_r10.fq.gz | samtools view -b | samtools sort -o $datadir/align_r10/r10.sorted.bam -T $datadir/align_r10/reads.tmp
    samtools index $datadir/align_r10/r10.sorted.bam
    
    minimap2 -a -x map-ont -t 36 $ref $datadir/fastqs/zymo_r9.fq.gz | samtools view -b | samtools sort -o $datadir/align_r9/r9.sorted.bam -T $datadir/align_r9/reads.tmp
    samtools index $datadir/align_r9/r9.sorted.bam
fi

if [ $1 == cons ] ; then
    ml samtools
    ml bcftools

    mkdir -p $datadir/cons
    
    ref=$datadir/Reference/ecoli_k12.fa
    
    ##bcftools mpileup -f $ref $datadir/align_r10/r10_ecoli_sub30.sorted.bam | bcftools call --ploidy 1 -mv -Ob -o $datadir/cons/r10_calls.bcf
    bcftools mpileup -f $ref $datadir/align_r10/r10_ecoli_sub10.sorted.bam | bcftools call --ploidy 1 -mv -Ob -o $datadir/cons/r10_sub10_calls.bcf 
    ##bcftools mpileup -f $ref $datadir/align_r9/r9_ecoli.sorted.bam | bcftools call --ploidy 1 -mv -Ob -o $datadir/cons/r9_calls.bcf 
    
fi

if [ $1 == bcfidx ] ; then
    ml bcftools
    ml samtools

    bcftools index $datadir/cons/r9_calls.bcf &
    bcftools index $datadir/cons/r10_calls.bcf &
    bcftools index $datadir/cons/r10_sub10_calls.bcf
fi

if [ $1 == seqcons ] ; then
    ml samtools

    bcftools view $datadir/cons/r10_calls.bcf | vcfutils.pl vcf2fq > $datadir/cons/r10_ecoli.cons.fq
    bcftools view $datadir/cons/r9_calls.bcf | vcfutils.pl vcf2fq > $datadir/cons/r9_ecoli.cons.fq

    seqtk seq -a $datadir/cons/r10_ecoli.cons.fq > $datadir/cons/r10_ecoli.cons.fa
    seqtk seq -a $datadir/cons/r9_ecoli.cons.fq > $datadir/cons/r9_ecoli.cons.fa
fi
