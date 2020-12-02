#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/nivar
datadir=/uru/Data/Nanopore/projects/nivar/paperfigs
fq=/uru/Data/Nanopore/projects/nivar/r9/r9_3kb.fq
dbxdir=~/Dropbox/yfan/nivar/paperfigs


ref=$rawdir/reference/candida_nivariensis.fa
gla=$rawdir/reference/medusa_fungi/candida_glabrata.fa

fin=$datadir/assembly_final/nivar.final.fasta
gff=$datadir/annotation_final/nivar.final.gff

if [ $1 == amino ] ; then
    ## ~/software/Augustus/scripts/gtf2aa.pl $fin $datadir/annotation/braker/braker.gtf $datadir/annotation_final/nivar.final.faa
     ~/software/Augustus/scripts/gtf2aa.pl $fin $gff $datadir/annotation_final/nivar.final.faa
fi

