#!/bin/bash

dbdir=/atium/Data/ref/mtases
gendir=$dbdir/genomes

if [ $1 == cleanstuff ] ; then
    ##from rebase ftp://ftp.neb.com/pub/rebase/Type_II_methyltransferase_genes_DNA.txt
    awk '{if (substr($1,1,1) ~ ">" ) print $1; else print $0}' $dbdir/Type_II_methyltransferase_genes_DNA.txt | sed -r 's/\s//g' > $dbdir/mtases.fa
fi

