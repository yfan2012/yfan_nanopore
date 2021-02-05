#!/bin/bash

dbdir=/mithril/Data/Nanopore/ref/kraken2

if [ $1 == builddb_ont ] ; then
    ##build db for ont
    kraken2-build --standard --threads 12 \
		  --kmer-len 28 \
		  --minimizer-len 24 \
		  --minimizer-spaces 6 \
		  --use-ftp \
		  --db $dbdir/standard_ont
fi

if [ $1 == builddb_illumina ] ; then
    ##build db for illumina
    kraken2-build --standard --threads 12 \
		  --use-ftp \
		  --db $dbdir/standard
fi

datadir=/uru/Data/Nanopore/projects/mdr/MDRstool_16

if [ $1 == runkraken_ont ] ; then
    ##run kraken2 on both 
    for i in native pcr ;
    do
	fq=$datadir/fastqs/MDRstool_16_$i.fq
	kraken2 \
	    --db $dbdir/standard_ont \
	    --threads 36 \
	    --classified-out $datadir/kraken/$i.class.txt \
	    --unclassified-out $datadir/kraken/$i.unclass.txt \
	    --output $datadir/kraken/$i.out.txt \
	    --report $datadir/kraken/$i.report.txt \
	    --use-names \
	    $fq
    done
fi

if [ $1 == runkraken_illumina ] ; then
    ##run kraken2 on illumina data
    for i in phase shotgun ;
    do
	fq=$datadir/illumina/trimmed/*$i*fq.gz
	kraken2 \
	    --db $dbdir/standard \
	    --threads 36 \
	    --classified-out $datadir/kraken/$i.class.txt \
	    --unclassified-out $datadir/kraken/$i.unclass.txt \
	    --output $datadir/kraken/$i.out.txt \
	    --report $datadir/kraken/$i.report.txt \
	    --use-names \
	    --gzip-compressed \
	    $fq
    done
fi
	    

if [ $1 == top40 ] ; then
    for i in phase shotgun pcr native ;
    do
	awk '$4 == "S" {print $0}' $datadir/kraken/$i.report.txt | \
	    sort -r -k2 -n | \
	    head -n 40 > $datadir/kraken/$i.report.top40.txt
    done
fi
    
    
