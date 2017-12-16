#!/bin/bash

penAdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs
ref=$penAdir/bcepacia_170816_BUCC.pilon.fa
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc

~/software/parsnp/parsnp -r $ref -d $penAdir -p 12 -o ${outdir}/bcc -c
