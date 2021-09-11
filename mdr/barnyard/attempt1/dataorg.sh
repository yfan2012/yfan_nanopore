#!/bin/bash

rawdir=/atium/Data/projects/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/barnyard
ssddir=~/data/mdr/barnyard

##prefix=210730_mdr_barnyard_mix1
##prefix=210730_mdr_barnyard_mix2
##prefix=210730_mdr_barnyard_mix3
prefix=210730_mdr_barnyard_mix4

if [ $1 == dlref ] ; then
    mkdir -p $datadir/ref
    ##plasmid from rachael
    ##ecoli
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O $datadir/ref/ecoli.fa.gz
    ##aureus
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz -O $datadir/ref/staph.fa.gz
fi

if [ $1 == catref ] ; then
    cat $datadir/ref/ecoli.fa.gz $datadir/ref/pJW792.fa.gz $datadir/ref/staph.fa.gz > $datadir/ref/allrefs.fa.gz
    gunzip $datadir/ref/allrefs.fa.gz
fi

if [ $1 == untar ] ; then
    mkdir -p $ssddir/raw
    tar -xzf $rawdir/$prefix.tar.gz -C $ssddir/raw/
fi

if [ $1 == slim ] ; then
    ##assumes tarball is appropriately exploded
    ##delete things that aren't in the subset so analysis doesn't take forever
    
    npid=`basename $ssddir/raw/$prefix/no_sample/*/fast5_pass/*_0.fast5 _0.fast5`
    for i in 0 1 2 3 4 ;
    do
	rm $ssddir/raw/$prefix/no_sample/*/fast5_pass/$npid*_$i*.fast5
    done

    for i in 5 6 7 8 9 ;
    do
	rm $ssddir/raw/$prefix/no_sample/*/fast5_pass/$npid*_$i.fast5
    done

    for i in {50..99} ;
    do
	rm $ssddir/raw/$prefix/no_sample/*/fast5_pass/$npid*_$i.fast5
    done
fi
	
if [ $1 == subset ] ; then
    ##grab subset of 2 million reads for toying with later, just in case it's needed
    ##do this after the 'slim' step above
    mkdir -p $datadir/raw_subset
    mkdir -p $datadir/raw_subset/$prefix

    cp $ssddir/raw/$prefix/no_sample/*/fast5_pass/*.fast5 $datadir/raw_subset/$prefix/
fi
   
