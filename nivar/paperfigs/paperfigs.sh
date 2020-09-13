#!/bin/bash

datadir=/uru/Data/Nanopore/projects/nivar

if [ $1 == medusa ] ; then
    scadir=$datadir/medusa
    mkdir -p $scadir

    cp -r ~/software/medusa/medusa_scripts ./
    
    java -jar ~/software/medusa/medusa.jar -f $datadir/reference/medusa_fungi -i $datadir/pilon/r9_pilon/nivar_r9.pilon_bwa.6.fasta -v -o $scadir/nivar_r9.pilon_bwa.6.scaffold.fasta

    java -jar ~/software/medusa/medusa.jar -f $datadir/reference/medusa_fungi -i $datadir/pilon/r10_pilon/nivar_r10.pilon_bwa.6.fasta -v -o $scadir/nivar_r10.pilon_bwa.6.scaffold.fasta
fi
