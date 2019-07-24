#!/bin/bash

datadir=~/work/seqlab/fungus

if [ $1 == rearrange ] ; then
    for run in run1 run2 ;
    do
	mkdir $datadir/$run/logs
	##move everything to different folders for batch calling
	for i in $datadir/$run/fast5_pass/*fast5 ;
	do
	    filename=`basename $i .fast5`
	    num=`echo $filename | cut -d _ -f 3`
	    mkdir -p $datadir/$run/fast5_pass/$num
	    mv $i $datadir/$run/fast5_pass/$num
	done
    done
fi

if [ $1 == submit ] ; then
    sbatch --array=0-145 --job-name=run1 --output=$datadir/run1/logs/%A_%a.out ./bc_call.scr $datadir/run1
    sbatch --array=0-247 --job-name=run2 --output=$datadir/run2/logs/%A_%a.out ./bc_call.scr $datadir/run2
fi


if [ $1 == cat_codes ] ; then
    for i in 01 02 03 04 05 06 07 08 09 10 11 12 ;
    do
	cat $datadir/run1/barcoded/*/barcode${i}/*fastq > $datadir/run1/barcode${i}.fastq
	cat $datadir/run2/barcoded/*/barcode${i}/*fastq > $datadir/run2/barcode${i}.fastq
    done
fi

if [ $1 == cat_fq ] ; then
    cat $datadir/run1/barcode09.fastq $datadir/run2/barcode09.fastq > $datadir/candida_nivariensis.fastq
fi


