#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/moth

if [ $1 == assemble ] ; then
   wtdbg2 -t 32 -i $datadir/181005_moth_cat.fastq -fo $datadir/wtdbg2
fi

if [ $1 == cons ] ; then
   wtpoa-cns -t 32 -i $datadir/wtdbg2.ctg.lay.gz -fo $datadir/wtdbg2.ctg.lay.fasta
fi
