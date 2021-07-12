#!/bin/bash

mtasedb=/atium/Data/ref/mtases

projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/plasmid

if [ $1 == makeblastdb ] ; then
    makeblastdb \
	-in $mtasedb/mtases.fa \
	-out $mtasedb/mtases \
	-dbtype nucl
fi

if [ $1 == blast ] ; then
    mkdir -p $datadir/blast

    for i in `cat rage_genomes.txt` ;
    do
	blastn \
	    -num_threads 36 \
	    -query $datadir/genomes/$i.fa \
	    -db $mtasedb/mtases \
	    -outfmt 7 \
	    -out $datadir/blast/$i.tsv
    done
fi

	
