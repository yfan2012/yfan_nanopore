#!/bin/bash

dbdir=/mithril/Data/Nanopore/ref
datadir=/dilithium/Data/Nanopore/projects/searslab_testont

prefix=1_7_million_psampleAPCflfl_Exp1

if [ $1 == stdbuild ] ; then
    
    bracken-build \
	-d $dbdir/kraken2/standard \
	-t 36 \
	-l 500
fi


if [ $1 == ontbuild ] ; then
    bracken-build \
	-d $dbdir/kraken2/standard_ont \
	-k 28 \
	-t 36 \
	-l 500
fi

    
     
