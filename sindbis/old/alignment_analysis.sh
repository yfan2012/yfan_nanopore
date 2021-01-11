#!/bin/bash

srcdir=~/Code/utils
outdir=~/Dropbox/timplab_data/sindbis
datadir=/dilithium/Data/Nanopore/sindbis

if [ $1 == genomecov ]; then
    ##coverage across the genome
    for i in mock infected antibody ;
    do
	mkdir -p $datadir/$i/cov
	bedtools genomecov -d -ibam $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.cov
    done
fi

if [ $1 == run_stats ] ; then
    ##get run stats in a nice csv, so you can grab this info for better normalization
    touch $outdir/runstats.csv
    echo file,yield,numreads,avglen >> $outdir/runstats.csv
    for i in mock infected antibody ;
    do
	bash $srcdir/qc/basic_run_assess.sh $datadir/$i/fqs/$i.fq >> $outdir/runstats.csv
    done
fi
if [ $1 == run_stats2 ] ; then
    ##get run stats in a nice csv, so you can grab this info for better normalization
    touch $outdir/runstats2.csv
    echo file,yield,numreads,avglen >> $outdir/runstats2.csv
    for i in Antibody_2dpi Antibody_3dpi Sindbis_2dpi Sindbis_3dpi ;
    do
	bash $srcdir/qc/basic_run_assess.sh $datadir/$i/fqs/$i.fq >> $outdir/runstats2.csv
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
