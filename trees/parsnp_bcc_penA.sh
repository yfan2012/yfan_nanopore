#!/bin/bash

penAdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_penA_genes
ref=$penAdir/penA_bcepacia_170816_BUCC_pilon.fasta
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc

~/software/parsnp/parsnp -r $ref -d $penAdir -p 12 -o ${outdir}/bcc_penA -c
