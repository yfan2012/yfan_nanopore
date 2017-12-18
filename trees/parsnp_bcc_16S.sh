#!/bin/bash

rrna16Sdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_16S_genes
ref=$rrna16Sdir/bcepacia_170816_BUCC_pilon0.fa
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc

~/software/parsnp/parsnp -r $ref -d $rrna16Sdir -p 12 -o ${outdir}/bcc_16S -c
