#!/bin/bash

rawdir=/atium/Data/projects/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard

prefix=210730_mdr_barnyard_mix4

if [ $1 == dlref ] ; then
    mkdir -p $datadir/ref
    
    ##already have ecoli
    ##already have staph
    
    ##strep
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/267/845/GCF_001267845.1_ASM126784v1/GCF_001267845.1_ASM126784v1_genomic.fna.gz -O $datadir/ref/strep.fa.gz
    
    ##bacilus
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz -O $datadir/ref/bacillus.fa.gz
fi

if [ $1 == catref ] ; then
    cat $datadir/ref/ecoli.fa.gz $datadir/ref/staph.fa.gz $datadir/ref/strep.fa.gz $datadir/ref/bacillus.fa.gz > $datadir/ref/speciescheck_refs.fa.gz
    gunzip $datadir/ref/speciescheck_refs.fa.gz
fi

   
