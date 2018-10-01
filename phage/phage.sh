#!/bin/bash

datadir=~/work/180809_phage
srcdir=~/Code/utils/marcc

if [ $1 == untar ] ; then
    mkdir -p $datadir/batch_logs
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_phage $srcdir/untar.scr $datadir/180809_bacteriophage_barcoded_15_20.raw.tgz $datadir
fi
    
if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done
    sbatch --array=0-149 --job-name=call_phage --output=$datadir/call_logs/phage.%A_%a.out $srcdir/bc_call.scr $datadir
fi
