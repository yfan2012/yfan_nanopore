#!/bin/bash

##datadir=/kyber/Data/seqlab/sp_2019/fungus_asm
datadir=~/work/seqlab/fungus

if [ $1 == maker_ctlfiles ] ; then
    maker -CTL
fi

if [ $1 == augustus ] ; then
    anndir=$datadir/annotation
    mkdir -p $anndir/augustus
    
    augustus --species=candida_albicans $anndir/candida_nivariensis_canu.pilon.6.fasta > $anndir/augustus/cani_by_caal.gff
fi

if [ $1 == medusa ] ; then
    scadir=$datadir/medusa
    mkdir -p $scadir

    ##cp $datadir/assembly/canu_pilon/candida_nivariensis_canu.pilon.6.fasta $scadir
    ##sed -i -e 's/_pilon//g' $scadir/candida_nivariensis_canu.pilon.6.fasta
    
    java -jar ~/scratch/software/medusa/medusa.jar -f ~/work/seqlab/fungus/References -i $scadir/candida_nivariensis_canu.pilon.6.fasta -o $scadir/candida_nivariensis_scaffold.fasta -v
fi
