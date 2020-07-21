#!/bin/bash

##use bacteria/virus/archea database
dbdir=~/scratch/centrifuge_db
srcdir=~/scratch/centrifuge
datadir=/scratch/groups/mschatz1/cpowgs/mdr5

if [ $1 == mdr5 ] ; then
    for i in mdrstool_16 mdrstool_19 recswab_1r ;
    do
	mkdir -p $datadir/$i/classification
	$srcdir/centrifuge -p 36 -x $dbdir/abv -U $datadir/$i/fastqs/$i.fq -S $datadir/$i/classification/$i.txt --report-file $datadir/$i/classification/${i}_report.tsv
	$srcdir/centrifuge-kreport -x $dbdir/abv $datadir/$i/classification/$i.txt > $datadir/$i/classification/kreport_$i.txt
    done
fi
