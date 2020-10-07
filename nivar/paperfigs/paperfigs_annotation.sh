#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/nivar
datadir=/uru/Data/Nanopore/projects/nivar/paperfigs
fq=/uru/Data/Nanopore/projects/nivar/r9/r9_3kb.fq
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
