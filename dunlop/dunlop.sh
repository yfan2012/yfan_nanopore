#!/bin/bash

datadir=~/work/180714_dunlop_3ecoli
srcdir=~/Code/utils/marcc

mkdir -p $datadir/batch_logs

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_dunlop $srcdir/untar.scr $datadir/180714_dunlop_3ecoli.tar.gz $datadir
fi

   
if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done
    sbatch --array=0-693 --job-name=call_dunlop --output=$datadir/call_logs/dunlop_ecoli3.%A_%a.out $srcdir/bc_call.scr $datadir
fi


if [ $1 == illumina_dilith ] ; then
    ##copy illumina data to dilithium
    mkdir -p /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli
    for i in `find ~/BaseSpace/Projects/180722_dunlop_3ecoli/ -name *fastq.gz` ;
    do
	cp $i /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli/
    done
fi

