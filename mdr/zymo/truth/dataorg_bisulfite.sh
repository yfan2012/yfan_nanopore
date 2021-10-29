#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/zymo/truth

if [ $1 == dl_wgbs ] ; then
    ##run info downloaded from ncbi web
    ##selected only bacterial wgbs runs
    mkdir -p $datadir/bisulfite
    awk -F ',' '{print $1}' $datadir/bisulfite/raw/wgbs_runinfo.csv > $datadir/bisulfite/raw/wgbs_acc.txt

    mkdir -p $datadir/bisulfite/raw/fastqs
    while read p ; do
	fastq-dump \
	    -O $datadir/bisulfite/raw/fastqs \
	    --split-files \
	    $p
    done < $datadir/bisulfite/raw/wgbs_acc.txt
fi

if [ $1 == rename ] ; then

    awk -F , '{print $1, $30}' $datadir/bisulfite/raw/wgbs_runinfo.csv > $datadir/bisulfite/raw/wgbs_names.txt

    while read p ; do
	acc=`echo $p | cut -d ' ' -f 1`
	name=`echo $p | cut -d ' ' -f 2`

	mv $datadir/bisulfite/raw/fastqs/${acc}_1.fastq $datadir/bisulfite/raw/fastqs/${name}_1.fastq
	mv $datadir/bisulfite/raw/fastqs/${acc}_2.fastq $datadir/bisulfite/raw/fastqs/${name}_2.fastq
    done < $datadir/bisulfite/raw/wgbs_names.txt
fi
