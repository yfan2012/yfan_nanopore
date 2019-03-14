#!/bin/bash

datadir=/kyber/Data/seqlab/sp_2019/fungus_asm

if [ $1 == maker_ctlfiles ] ; then
    maker -CTL
fi

if [ $1 == augustus ] ; then
    anndir=$datadir/annotation
    mkdir -p $anndir/augustus
    
    augustus --species=candida_albicans $anndir/candida_nivariensis_canu.pilon.6.fasta > $anndir/augustus/cani_by_caal.gff
fi

