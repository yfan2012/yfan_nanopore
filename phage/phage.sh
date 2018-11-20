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

if [ $1 == fastqs ] ; then
    mkdir -p $datadir/fastqs
    cat $datadir/called/*/workspace/pass/barcode01/*.fastq > $datadir/fastqs/samp15at24hrs.fq
    cat $datadir/called/*/workspace/pass/barcode02/*.fastq > $datadir/fastqs/samp15at48hrs.fq
    cat $datadir/called/*/workspace/pass/barcode03/*.fastq > $datadir/fastqs/samp20at24hrs.fq
    cat $datadir/called/*/workspace/pass/barcode04/*.fastq > $datadir/fastqs/samp20at48hrs.fq
fi


if [ $1 == vanilla_centrifuge ] ; then
    mkdir -p $datadir/centrifuge
    
