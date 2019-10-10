#!/bin/bash

##script for everything up to assembly
datadir=/uru/Data/Nanopore/projects/nivar

if [ $1 == r10_unpack ] ; then
    prefix=r10
    
    tarball=/kyber/Data/Nanopore/oxford/190706_nivar_r10/190706_nivar_r10.tar.gz

    rawdir=$datadir/$prefix/raw
    mkdir -p $rawdir
    
    tar -xzf $tarball -C $datadir/$prefix/raw
fi


if [ $1 == r9_call ] ; then
    raw=/uru/Data/Nanopore/projects/nivar/r9/raw
    guppy_basecaller -i $raw/run1 -o $raw/run1_called -c dna_r9.4.1_450bps -x "cuda:0"
    guppy_basecaller -i $raw/run2 -o $raw/run2_called -c dna_r9.4.1_450bps -x "cuda:0"

fi


