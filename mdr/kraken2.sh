#!/bin/bash

dbdir=/mithril/Data/Nanopore/ref/kraken2

if [ $1 == builddb ] ; then
    kraken2-build --standard --threads 12 --kmer-len 28 --minimizer-len 24 --minimizer-spaces 6 --use-ftp --db $dbdir/standard
    ##https://github.com/DerrickWood/kraken2/issues/226
fi

datadir=/uru/Data/Nanopore/projects/mdr/MDRstool_16

if [ $1 == runkraken ] ; then
    natfq=$datadir/fastqs/MDRstool_16_native.fq
    pcrfq=$datadir/fastqs/MDRstool_16_pcr.fq
    illfq=$datadir
	
    
