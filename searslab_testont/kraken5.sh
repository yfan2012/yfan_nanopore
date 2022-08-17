#!/bin/bash

datadir=/dilithium/Data/Nanopore/projects/searslab_testont
dbdir=/mithril/Data/Nanopore/ref/kraken2
db=$dbdir/standard
ontdb=$dbdir/standard_ont

run=5_millions_reads_psample_APCflfl_Exp1

if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs
    mkdir -p $datadir/fastqs/$run

    for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ;
    do
	cat $datadir/raw/$run/no_sample/*/fastq_pass/barcode$i/*fastq.gz \
	    > $datadir/fastqs/$run/barcode$i.fastq.gz
    done
fi

       
if [ $1 == ontkraken ] ; then
    mkdir -p $datadir/ontkraken
    mkdir -p $datadir/ontkraken/$run

    for bc in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ;
    do
	i=barcode$bc
	kraken2 \
	    --db $ontdb \
	    --threads 12 \
	    --classified-out $datadir/kraken/$run/$i.class.txt \
	    --unclassified-out $datadir/kraken/$run/$i.unclass.txt \
	    --output $datadir/kraken/$run/$i.out.txt \
	    --report $datadir/kraken/$run/$i.report.txt \
	    --use-names \
	    --gzip-compressed \
	    $datadir/fastqs/$run/$i.fastq.gz
    done
fi


if [ $1 == stdkraken ] ; then
    mkdir -p $datadir/stdkraken
    mkdir -p $datadir/stdkraken/$run

    for bc in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ;
    do
	i=barcode$bc
	kraken2 \
	    --db $db \
	    --threads 12 \
	    --classified-out $datadir/stdkraken/$run/$i.class.txt \
	    --unclassified-out $datadir/stdkraken/$run/$i.unclass.txt \
	    --output $datadir/stdkraken/$run/$i.out.txt \
	    --report $datadir/stdkraken/$run/$i.report.txt \
	    --use-names \
	    --gzip-compressed \
	    $datadir/fastqs/$run/$i.fastq.gz
    done
fi
