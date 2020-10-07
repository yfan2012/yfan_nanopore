#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/nivar
datadir=$rawdir/paperfigs
fq=$rawdir/r9/r9_3kb.fq
drna=$rawdir/dRNA/dRNA.fq
rnafwd=$rawdir/illumina/cDNA_trimmed/nivar_cDNA_fwd_paired.fq
rnarev=$rawdir/illumina/cDNA_trimmed/nivar_cDNA_rev_paired.fq
dbxdir=~/Dropbox/yfan/nivar/paperfigs

gen=$datadir/assembly_final/nivar.final.fasta
ref=$rawdir/reference/candida_nivariensis.fa

gla=$rawdir/reference/medusa_fungi/candida_glabrata.fa
cer=$rawdir/reference/saccharomyces_cerevisiae.fa
cergff=$rawdir/reference/saccharomyces_cerevisiae.gff

if [ $1 == liftoff_cerevisiae ] ; then
    mkdir -p $datadir/annotation/liftoff
    liftoff \
	-p 36 \
	-o $datadir/annotation/liftoff/nivar_cer_lifted.gff \
	-u $datadir/annotation/liftoff/nivar_cer_unmapped.txt \
	-dir $datadir/annotation/liftoff/intermediate_files \
	-g $cergff \
	$gen \
	$cer
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/annotation/align

    minimap2 -t 36 \
	-ax splice -uf -k14 $gen $drna | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/annotation/align/nivar_dRNA.sorted.bam
    samtools index $datadir/annotation/align/nivar_dRNA.sorted.bam

    echo hisat####################################
    ##hisat2-build $gen $datadir/assembly_final/nivar.final
    hisat2 -p 36 \
	-x $datadir/assembly_final/nivar.final \
	-1 $rnafwd \
	-2 $rnarev | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/annotation/align/nivar_cDNA.sorted.bam
    samtools index $datadir/annotation/align/nivar_cDNA.sorted.bam
fi

if [ $1 == intersect ] ; then
    ##check how many liftoff genes have support
