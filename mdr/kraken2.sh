#!/bin/bash

dbdir=/mithril/Data/Nanopore/ref/kraken2

if [ $1 == builddb ] ; then
    kraken2-build --standard --threads 12 --db $dbdir/standard
    ##https://github.com/DerrickWood/kraken2/issues/226
fi
