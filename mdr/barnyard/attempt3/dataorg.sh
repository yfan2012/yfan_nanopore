#!/bin/bash

rawdir=/atium/Data/projects/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/barnyard

if [ $1 == dlref ] ; then
    mkdir -p $datadir/ref
    ##plasmid from rachael
    ##ecoli
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O $datadir/ref/ecoli.fa.gz
    ##aureus
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz -O $datadir/ref/staph.fa.gz
fi

if [ $1 == catref ] ; then
    cat $datadir/ref/dhq5a.fa.gz $datadir/ref/pRW62.fa.gz $datadir/ref/staph.fa.gz > $datadir/ref/allrefs.fa.gz
    gunzip $datadir/ref/allrefs.fa.gz
fi

   
