#!/bin/bash

datadir=/uru/Data/Nanopore/projects/nivar
dbxdir=~/Dropbox/yfan/nivar

if [ $1 == find_enrich ] ; then
    
    ref=$datadir/reference/candida_nivariensis.fa
    mkdir -p $dbxdir/motif_enrich
    
    for i in $datadir/mummer/*ref/*.snps ;
    do
	prefix=`basename $i .snps`
	echo $prefix
	python2 ~/Code/utils/motif_enrich.py -s $i -r $ref -m 6 -o $dbxdir/motif_enrich/$prefix.csv
    done

    for i in $datadir/mummer/*raw/*.snps ;
    do
	prefix=`basename $i .snps`
	echo $prefix
	corr=`echo $prefix | cut -d _ -f 3 `
	pore=`echo $prefix | cut -d _ -f 2`

	refcorr=$datadir/mummer/${pore}_${corr}_raw/$corr.15.fasta

	python2 ~/Code/utils/motif_enrich.py -s $i -r $refcorr -m 6 -o $dbxdir/motif_enrich/$prefix.csv
    done

fi




    
