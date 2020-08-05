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
	prefix=`echo $fq | basename .fq`
	kraken2 \
	    --db $dbdir/standard_ont \
	    --threads 36 \
	    --classified-out $datadir/kraken/$prefix.class.txt \
	    --unclassified-out $datadir/kraken/$prefix.unclass.txt \
	    --output $datadir/kraken/$prefix.out.txt \
	    --report $datadir/kraken/$prefix.report.txt \
	    --use-names \
	    $fq
    done
fi

if [ $1 == runkraken_illumina ] ; then
    ##run kraken2 on illumina data
    for i in phase shotgun ;
    do
	fq=$datadir/illumina/*$i*fq.gz
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
	    

    
    
    
