#!/bin/bash

datadir=/uru/Data/Nanopore/projects/nivar

ref=$datadir/reference/candida_nivariensis.fa
gla=$datadir/reference/candida_glabrata.fa
cer=$datadir/reference/saccharomyces_cerevisiae.fa
alb=$datadir/reference/candida_albicans.fa


fin=$datadir/paperfigs/assembly_final/nivar.final.fasta

fqraw=$datadir/r9/r9.fq
fq=$datadir/r9/r9_3kb.fq
illfwd=$datadir/illumina/gDNA_trimmed/nivar_gDNA_fwd_paired.fq.gz
illrev=$datadir/illumina/gDNA_trimmed/nivar_gDNA_rev_paired.fq.gz

if [ $1 == trf ] ; then

    mkdir -p $datadir/paperfigs/trf
    ##File Match Mismatch Delta PM PI Minscore MaxPeriod [options]
    ~/software/trf $fin 2 7 7 80 10 50 600 -ngs
fi

    
