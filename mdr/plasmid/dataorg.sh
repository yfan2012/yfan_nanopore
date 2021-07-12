#!/bin/bash

mtasedb=/atium/Data/ref/mtases
projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/plasmid

if [ $1 == dlgenomes ] ; then
    mkdir -p $datadir/genomes/$p.fa
    for p in `cat rage_genomes.txt` ;
    do
	echo $p
	esearch -db nuccore -query $p | efetch -format fasta > $datadir/genomes/$p.fa
    done 
    
fi
