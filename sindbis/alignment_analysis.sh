#!/bin/bash

##datadir=/scratch/groups/mschatz1/cpowgs/sindbis
datadir=~/Dropbox/Timplab_Data/sindbis

if [ $1 == genomecov ]; then
    ##coverage across the genome
    for i in mock infected antibody ;
    do
	mkdir -p $datadir/$i/cov
	bedtools genomecov -d -ibam $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.cov
    done
fi

if [ $1 == genomecov_v2 ]; then
    ##coverage across the genome
    for i in Antibody_2dpi Antibody_3dpi Sindbis_2dpi Sindbis_3dpi ;
    do
	mkdir -p $datadir/$i/cov
	bedtools genomecov -d -ibam $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.cov &
    done
fi
