#!/bin/bash

penAdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/penA_genes
ref=$penAdir/penA_170816_BUCC.pilon.fasta
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc

~/software/parsnp/parsnp -r $ref -d $penAdir -p 12 -o ${outdir}/penA -c
