#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/sindbis

if [ $1 == genomecov ]; then
    ##coverage across the genome
    for i in mock infected antibody ;
    do
	mkdir -p $datadir/$i/cov
	bedtools genomecov -d -ibam $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.cov
    done
fi
