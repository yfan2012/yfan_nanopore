#!/bin/bash

rawdir=/work-zfs/mschatz1/cpowgs
datadir=/scratch/groups/mschatz1/cpowgs/BOAV
srcdir=~/Code/utils/marcc

mkdir -p $datadir

if [ $1 == untar ] ; then
    mkdir -p $datadir/batch_logs
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar_log.txt $srcdir/untar.scr $rawdir/BOAV.tar.gz $datadir
fi

if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs

    numdirs=`find $datadir/raw/* -maxdepth 0 -type d | wc -l `
    dummy=1
    maxdir=`expr $numdirs - $dummy`

    sbatch --array=0-$maxdir --job-name=BOAV_call --output=$datadir/call_logs/BOAV.%A_%a.out $srcdir/call.scr $datadir
fi


if [ $1 == fq ] ; then
    mkdir -p $datadir/fastqs

    sbatch --output=$datadir/batch_logs/fq.out $srcdir/fqs.scr $datadir
fi


if [ $1 == assemble ] ; then
    mkdir -p $datadir/canu_assembly

    ##canu complaining about read qualities - copied assembly script and altered to not stop on read qual
    bash ./assemble_bacteria.sh $datadir/fastqs/BOAV.fq $datadir/canu_assembly
fi


if [ $1 == assemble17 ] ; then
    mkdir -p $datadir/canu17

    ##canu complaining about read qualities - copied assembly script and altered to not stop on read qual
    bash ./assemble_bacteria17.sh $datadir/fastqs/BOAV.fq $datadir/canu17
fi

if [ $1 == pilon ] ; then
    mkdir -p $datadir/pilon

    cp $rawdir/B-avium_S5_L001_R* $datadir/pilon
    rename B-avium BOAV $datadir/pilon/*
    sbatch --output=$datadir/batch_logs/pilon.out --job-name=BOAV $srcdir/pilon.scr $datadir
fi


if [ $1 == centrifuge ] ; then
    centdir=~/scratch/centrifuge
    dbdir=~/scratch/centrifuge_db

    mkdir -p $datadir/classification
    
    $centdir/centrifuge -p 36 -x $dbdir/abv -U $datadir/fastqs/BOAV.fq -S $datadir/classification/BOAV.txt --report-file $datadir/classification/BOAV_report.tsv
    $centdir/centrifuge-kreport -x $dbdir/abv $datadir/classification/BOAV.txt > $datadir/classification/kreport_BOAV.txt
fi

if [ $1 == illumina_centrifuge ] ; then
    centdir=~/scratch/centrifuge
    dbdir=~/scratch/centrifuge_db

    mkdir -p $datadir/classification_illumina
    
    $centdir/centrifuge -p 36 -x $dbdir/abv -1 $rawdir/B-avium_S5_L001_R1_001.fastq.gz -2 $rawdir/B-avium_S5_L001_R2_001.fastq.gz -S $datadir/classification/illumina_BOAV.txt --report-file $datadir/classification/illumina_BOAV_report.tsv
    $centdir/centrifuge-kreport -x $dbdir/abv $datadir/classification/illumina_BOAV.txt > $datadir/classification/kreport_illumina_BOAV.txt
fi

