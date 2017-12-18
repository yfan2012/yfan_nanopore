#!/bin/bash

datadir=/atium/Data/Nanopore/cpowgs/170816_BUCC
refdir=$datadir/bcc_refs
rrna16Sdir=$datadir/bcc_16S

for i in $refdir/*fa ;
do
    fafile=`echo $i | cut -d '/' -f 8`
    prefix=`echo $fafile | cut -d '.' -f 1`
    barrnap --quiet --threads 12 $i > $rrna16Sdir/$prefix.gff
    python ~/Code/carbapenem_r21/bucc/16S_extract.py -g $rrna16Sdir/$prefix.gff
done

    
