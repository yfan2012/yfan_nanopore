#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/mdr5/CAAU_1

if [ $1 == assemble ] ; then
    mkdir -p $datadir/wtdbg2
    wtdbg2 -t 32 -p 19 -A -S 1 -s 0.05 -K 1000 -L 5000 -i $datadir/fastqs/CAAU_1.fq -fo $datadir/wtdbg2/CAAU_1
fi

if [ $1 == cons ] ; then
    wtpoa-cns -t 32 -i $datadir/wtdbg2/CAAU_1.ctg.lay.gz -fo $datadir/wtdbg2/CAAU_1.ctg.lay.fasta
fi


