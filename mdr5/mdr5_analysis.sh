#!/bin/bash

datadir=/dilithium/Data/Nanopore/mdr5
abriout=~/Dropbox/Timplab_Data/mdr5/abricate


if [ $1 == abricate ] ; then
    ##abricate with card only for now
    mkdir -p $abriout
    for i in $datadir/assemblies/raw/*.fasta ;
    do
	prefix=`basename $i .contigs.fasta`
	abricate --db card $i > $abriout/$prefix.card.tsv
    done
fi

if [ $1 == abricate_report ] ; then
    for i in $abriout/*.card.tsv ;
    do
	prefix=`basename $i .card.tsv`
	python ~/Code/carbapenem_r21/report/abricate_report.py -p $abriout/$prefix -o $abriout/$prefix.abricate.csv
    done
fi
