#!/bin/bash

refdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc/wga

for i in $refdir/*.fa ;
do
    prefix=` echo ${i%.fa} | cut -d '/' -f 8 `

    ##nucmer -p $outdir/$prefix $refdir/bcontaminans.fa $i
    show-coords -T $outdir/$prefix.delta > $outdir/$prefix.tsv 
    
done
